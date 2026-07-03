//! VcfStore — one-shot VCF preprocessing + region query.
//!
//! (docs/design.md §5.2) The VCF is the heaviest raw material. We parse it
//! once at load time (eager) into an in-memory representation that does NOT
//! assume biallelic / fixed ploidy, then serve region queries from memory.
//!
//! GT decoding uses `bcf_get_format_values` (the impl behind the
//! `bcf_get_genotypes` macro) so htslib normalises FORMAT GT width
//! (INT8/INT16/INT32) into a flat encoded int32 array; we only replicate the
//! `bcf_gt_*` bit operations in Rust (see `htslib_ffi`).

use crate::compiled::{
    AlleleOp, AlleleOpKind, CompiledRecord, RecordFlags, RecordKind, VcfCompileStats,
};
use crate::htslib_ffi as ffi;
use crate::planner::{plan_region, plan_region_set, PlanOptions, RegionPlan};
use smallvec::SmallVec;
use std::collections::HashMap;
use std::ffi::CString;
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::os::raw::{c_int, c_void};
use std::path::{Path, PathBuf};
use std::time::UNIX_EPOCH;

const CVCF_MAGIC: &[u8; 8] = b"CVCF0001";
const CVCF_VERSION: u32 = 6;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct SourceFingerprint {
    len: u64,
    mtime_secs: i64,
    mtime_nanos: u32,
}

/// One allele of a sample's genotype.
#[derive(Clone, Debug)]
pub struct GtAllele {
    /// `None` = missing (`.`); `Some(0)` = REF, `Some(1..)` = ALT index.
    pub allele: Option<i32>,
    /// `bcf_gt_is_phased(raw)` — phase bit of this allele.
    pub phased: bool,
    /// Original htslib-encoded GT int32, kept for exact diffing / diagnostics.
    pub raw: i32,
}

const COMPACT_GT_ALLELE_MASK: u16 = 0x3fff;
const COMPACT_GT_MISSING: u16 = 0x4000;
const COMPACT_GT_PHASED: u16 = 0x8000;

/// Compact general GT store for non-bitset selection paths.
///
/// Each GT allele is encoded in one u16: lower 14 bits store allele index,
/// bit 14 marks missing, and bit 15 stores the phased flag. Normal haplotype
/// fallback selection walks this contiguous representation; the legacy
/// `GtAllele` matrix is kept only when an allele index cannot fit here.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CompactGt {
    sample_offsets: Vec<u32>,
    codes: Vec<u16>,
}

impl CompactGt {
    pub fn from_gt(n_samples: usize, gt: &[SmallVec<[GtAllele; 2]>]) -> Option<Self> {
        if n_samples == 0 || gt.is_empty() || gt.len() != n_samples {
            return None;
        }
        let mut sample_offsets = Vec::with_capacity(n_samples + 1);
        let mut codes = Vec::new();
        sample_offsets.push(0);
        for sample_gt in gt {
            for allele in sample_gt {
                let mut code = if let Some(idx) = allele.allele {
                    if idx < 0 || idx > COMPACT_GT_ALLELE_MASK as i32 {
                        return None;
                    }
                    idx as u16
                } else {
                    COMPACT_GT_MISSING
                };
                if allele.phased {
                    code |= COMPACT_GT_PHASED;
                }
                codes.push(code);
            }
            let offset = u32::try_from(codes.len()).ok()?;
            sample_offsets.push(offset);
        }
        Some(CompactGt {
            sample_offsets,
            codes,
        })
    }

    #[inline]
    pub fn n_samples(&self) -> usize {
        self.sample_offsets.len().saturating_sub(1)
    }

    #[inline]
    pub fn sample(&self, sample_idx: usize) -> Option<&[u16]> {
        let start = *self.sample_offsets.get(sample_idx)? as usize;
        let end = *self.sample_offsets.get(sample_idx + 1)? as usize;
        self.codes.get(start..end)
    }

    #[inline]
    pub fn allele(code: u16) -> Option<i32> {
        if code & COMPACT_GT_MISSING != 0 {
            None
        } else {
            Some((code & COMPACT_GT_ALLELE_MASK) as i32)
        }
    }

    #[inline]
    pub fn phased(code: u16) -> bool {
        code & COMPACT_GT_PHASED != 0
    }
}

/// Fast genotype representation for biallelic diploid records.
///
/// This is intentionally record-local for now: it gives the hot `-H 1/2`
/// selection path a branch-light bit test while preserving the existing general
/// GT structure for fallback and parity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BiallelicPhasedGtBits {
    n_samples: usize,
    hap1_alt: Vec<u64>,
    hap2_alt: Vec<u64>,
    missing: Vec<u64>,
}

impl BiallelicPhasedGtBits {
    pub fn from_gt(
        n_allele: usize,
        n_samples: usize,
        gt: &[SmallVec<[GtAllele; 2]>],
    ) -> Option<Self> {
        if n_allele != 2 || n_samples == 0 || gt.len() != n_samples {
            return None;
        }
        let n_words = n_samples.div_ceil(64);
        let mut bits = BiallelicPhasedGtBits {
            n_samples,
            hap1_alt: vec![0; n_words],
            hap2_alt: vec![0; n_words],
            missing: vec![0; n_words],
        };
        for (sample, sample_gt) in gt.iter().enumerate() {
            if sample_gt.len() != 2 {
                return None;
            }
            let word = sample / 64;
            let bit = 1u64 << (sample & 63);
            let a0 = &sample_gt[0];
            let a1 = &sample_gt[1];
            let Some(h1) = a0.allele else {
                bits.missing[word] |= bit;
                continue;
            };
            let Some(h2) = a1.allele else {
                bits.missing[word] |= bit;
                continue;
            };
            if !a0.phased || !a1.phased || !(0..=1).contains(&h1) || !(0..=1).contains(&h2) {
                return None;
            }
            if h1 == 1 {
                bits.hap1_alt[word] |= bit;
            }
            if h2 == 1 {
                bits.hap2_alt[word] |= bit;
            }
        }
        Some(bits)
    }

    #[inline]
    pub fn allele_for_hap(&self, sample_idx: usize, hap: usize) -> Option<Option<i32>> {
        if sample_idx >= self.n_samples || hap == 0 || hap > 2 {
            return None;
        }
        let word = sample_idx / 64;
        let bit = 1u64 << (sample_idx & 63);
        if self.missing[word] & bit != 0 {
            return Some(None);
        }
        let alt_bits = if hap == 1 {
            &self.hap1_alt
        } else {
            &self.hap2_alt
        };
        Some(Some(if alt_bits[word] & bit != 0 { 1 } else { 0 }))
    }

    #[inline]
    pub fn is_alt_for_hap(&self, sample_idx: usize, hap: usize) -> Option<bool> {
        match self.allele_for_hap(sample_idx, hap)? {
            Some(0) => Some(false),
            Some(1) => Some(true),
            _ => None,
        }
    }

    #[inline]
    pub fn alt_words_for_hap(&self, hap: usize) -> Option<&[u64]> {
        match hap {
            1 => Some(&self.hap1_alt),
            2 => Some(&self.hap2_alt),
            _ => None,
        }
    }

    #[inline]
    pub fn missing_words(&self) -> &[u64] {
        &self.missing
    }
}

/// A preprocessed VCF record. Alleles are REF + ALTs; GT is per-sample.
#[derive(Clone, Debug)]
pub struct VcfRecord {
    /// 0-based POS (`bcf1_t.pos`).
    pub pos: i64,
    /// REF length (`bcf1_t.rlen`).
    pub rlen: i32,
    /// Contig id within this VCF's header.
    pub rid: i32,
    /// REF + ALTs, variable count, each variable length.
    pub alleles: Vec<SmallVec<[u8; 16]>>,
    /// Rare raw-GT fallback. Common records use `gt_compact`/`gt_bits` and keep
    /// this empty to avoid carrying a cache-cold matrix through runtime.
    pub gt: Vec<SmallVec<[GtAllele; 2]>>,
    /// Compact GT fallback store for general selection modes.
    pub gt_compact: Option<CompactGt>,
    /// Bitset GT fastpath for biallelic phased diploid records.
    pub gt_bits: Option<BiallelicPhasedGtBits>,
    /// `bcf_get_variant_types` bitmask (VCF_SNP|MNP|INDEL|...), precomputed.
    pub var_type: i32,
    /// Preclassified record/allele metadata used by fastpath dispatch.
    pub compiled: CompiledRecord,
}

pub enum RecordSet<'a> {
    RefSlice(&'a [&'a VcfRecord]),
    IndexSlice {
        records: &'a [VcfRecord],
        idx: &'a [u32],
    },
    IndexFilteredPrefixAndSlice {
        records: &'a [VcfRecord],
        prefix_idx: &'a [u32],
        idx: &'a [u32],
        start: i64,
        end: i64,
        prefix_len: usize,
    },
    Empty,
}

pub enum RecordSetIter<'a, 's> {
    RefSlice(std::iter::Copied<std::slice::Iter<'s, &'a VcfRecord>>),
    IndexSlice {
        records: &'a [VcfRecord],
        iter: std::slice::Iter<'s, u32>,
    },
    IndexFilteredPrefixAndSlice {
        records: &'a [VcfRecord],
        prefix_iter: std::slice::Iter<'s, u32>,
        idx_iter: std::slice::Iter<'s, u32>,
        start: i64,
        end: i64,
        prefix_remaining: usize,
    },
    Empty,
}

impl<'a> RecordSet<'a> {
    pub fn from_ref_slice(records: &'a [&'a VcfRecord]) -> Self {
        if records.is_empty() {
            RecordSet::Empty
        } else {
            RecordSet::RefSlice(records)
        }
    }

    #[inline]
    pub fn len(&self) -> usize {
        match self {
            RecordSet::RefSlice(records) => records.len(),
            RecordSet::IndexSlice { idx, .. } => idx.len(),
            RecordSet::IndexFilteredPrefixAndSlice {
                prefix_len, idx, ..
            } => prefix_len + idx.len(),
            RecordSet::Empty => 0,
        }
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn iter(&self) -> RecordSetIter<'a, '_> {
        match self {
            RecordSet::RefSlice(records) => RecordSetIter::RefSlice(records.iter().copied()),
            RecordSet::IndexSlice { records, idx } => RecordSetIter::IndexSlice {
                records,
                iter: idx.iter(),
            },
            RecordSet::IndexFilteredPrefixAndSlice {
                records,
                prefix_idx,
                idx,
                start,
                end,
                prefix_len,
            } => RecordSetIter::IndexFilteredPrefixAndSlice {
                records,
                prefix_iter: prefix_idx.iter(),
                idx_iter: idx.iter(),
                start: *start,
                end: *end,
                prefix_remaining: *prefix_len,
            },
            RecordSet::Empty => RecordSetIter::Empty,
        }
    }

    pub fn to_refs(&self) -> Vec<&'a VcfRecord> {
        self.iter().collect()
    }
}

impl<'a> Iterator for RecordSetIter<'a, '_> {
    type Item = &'a VcfRecord;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            RecordSetIter::RefSlice(iter) => iter.next(),
            RecordSetIter::IndexSlice { records, iter } => {
                iter.next().map(|&i| &records[i as usize])
            }
            RecordSetIter::IndexFilteredPrefixAndSlice {
                records,
                prefix_iter,
                idx_iter,
                start,
                end,
                prefix_remaining,
            } => {
                for &i in prefix_iter.by_ref() {
                    let rec = &records[i as usize];
                    if rec.pos <= *end && rec.ref_end() >= *start {
                        *prefix_remaining -= 1;
                        return Some(rec);
                    }
                }
                *prefix_remaining = 0;
                idx_iter.next().map(|&i| &records[i as usize])
            }
            RecordSetIter::Empty => None,
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        match self {
            RecordSetIter::RefSlice(iter) => iter.size_hint(),
            RecordSetIter::IndexSlice { iter, .. } => iter.size_hint(),
            RecordSetIter::IndexFilteredPrefixAndSlice {
                idx_iter,
                prefix_remaining,
                ..
            } => {
                let (i_lo, i_hi) = idx_iter.size_hint();
                (
                    *prefix_remaining + i_lo,
                    i_hi.map(|i| *prefix_remaining + i),
                )
            }
            RecordSetIter::Empty => (0, Some(0)),
        }
    }
}

impl ExactSizeIterator for RecordSetIter<'_, '_> {}

impl VcfRecord {
    /// 0-based inclusive end of the REF span: `pos + rlen - 1`.
    pub fn ref_end(&self) -> i64 {
        self.pos + self.rlen as i64 - 1
    }

    pub fn is_snp(&self) -> bool {
        self.var_type & ffi::VCF_SNP != 0
    }
    pub fn is_indel(&self) -> bool {
        self.var_type & ffi::VCF_INDEL != 0
    }
    pub fn is_mnp(&self) -> bool {
        self.var_type & ffi::VCF_MNP != 0
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LoadStrategy {
    /// Preprocess the whole VCF into memory up front.
    Eager,
    /// Preprocess on first use of this VCF (M5 wires the lazy trigger).
    Lazy,
}

/// In-memory store of one VCF's preprocessed records + region index.
pub struct VcfStore {
    path: PathBuf,
    records: Vec<VcfRecord>,
    /// rid -> record indices sorted by pos (stable on file order).
    by_rid: HashMap<i32, Vec<u32>>,
    /// rid -> prefix-max of `ref_end` aligned with `by_rid`, so a query can
    /// skip the leading prefix that cannot reach back into the region.
    pmax_end: HashMap<i32, Vec<i64>>,
    /// contig name -> rid
    seq_names: HashMap<String, i32>,
    sample_names: Vec<String>,
    /// name -> sample index
    sample_idx: HashMap<String, i32>,
    has_gt: bool,
    n_sample: i32,
    compile_stats: VcfCompileStats,
}

impl VcfStore {
    /// Eagerly load + preprocess a VCF (BCF or VCF, plain or bgzipped).
    pub fn load(path: impl Into<PathBuf>) -> Result<Self, String> {
        Self::load_with_strategy(path, LoadStrategy::Eager)
    }

    pub fn load_with_strategy(
        path: impl Into<PathBuf>,
        strat: LoadStrategy,
    ) -> Result<Self, String> {
        let path = path.into();
        // LAZY is handled at the engine level (M5); here we always parse when
        // called. The flag is accepted for API completeness.
        let _ = strat;
        if let Some(store) = Self::try_load_default_cache(&path) {
            return Ok(store);
        }

        let mut store = VcfStore::empty(path.clone());
        store.parse()?;
        let _ = store.write_default_cache();
        Ok(store)
    }

    fn empty(path: PathBuf) -> Self {
        VcfStore {
            path,
            records: Vec::new(),
            by_rid: HashMap::new(),
            pmax_end: HashMap::new(),
            seq_names: HashMap::new(),
            sample_names: Vec::new(),
            sample_idx: HashMap::new(),
            has_gt: false,
            n_sample: 0,
            compile_stats: VcfCompileStats::default(),
        }
    }

    pub fn path(&self) -> &std::path::Path {
        &self.path
    }

    pub fn n_records(&self) -> usize {
        self.records.len()
    }

    pub fn n_sample(&self) -> i32 {
        self.n_sample
    }

    pub fn has_gt(&self) -> bool {
        self.has_gt
    }

    pub fn sample_names(&self) -> &[String] {
        &self.sample_names
    }

    pub fn sample_index(&self, name: &str) -> Option<i32> {
        self.sample_idx.get(name).copied()
    }

    pub fn rid_of(&self, chr: &str) -> Option<i32> {
        self.seq_names.get(chr).copied()
    }

    pub fn records(&self) -> &[VcfRecord] {
        &self.records
    }

    pub fn compile_stats(&self) -> &VcfCompileStats {
        &self.compile_stats
    }

    /// Query records overlapping `[start, end]` (0-based, inclusive), sorted by
    /// pos. `overlap`: 0 = POS in region; 1 = record span overlaps; 2 = variant
    /// span overlaps (approximated as record span for biallelic SNP/indel; see
    /// docs §3). Never misses a deletion/MNP whose `pos < start` but spans in.
    pub fn query(&self, chr: &str, start: i64, end: i64, overlap: u8) -> Vec<&VcfRecord> {
        self.query_set(chr, start, end, overlap).to_refs()
    }

    pub fn query_set(&self, chr: &str, start: i64, end: i64, overlap: u8) -> RecordSet<'_> {
        let rid = match self.seq_names.get(chr) {
            Some(r) => *r,
            None => return RecordSet::Empty,
        };
        let idx = match self.by_rid.get(&rid) {
            Some(v) => v,
            None => return RecordSet::Empty,
        };
        if idx.is_empty() {
            return RecordSet::Empty;
        }
        let pmax = self
            .pmax_end
            .get(&rid)
            .expect("pmax_end aligned with by_rid");

        // hi = first index whose record.pos > end
        let hi = idx.partition_point(|&i| self.records[i as usize].pos <= end);
        // lo_pos = first index whose record.pos >= start
        let lo_pos = idx.partition_point(|&i| self.records[i as usize].pos < start);

        let first_spanning = if (overlap == 1 || overlap == 2) && lo_pos > 0 {
            pmax[..lo_pos].partition_point(|&m| m < start)
        } else {
            lo_pos
        };
        match overlap {
            0 => {
                if lo_pos == hi {
                    RecordSet::Empty
                } else {
                    RecordSet::IndexSlice {
                        records: &self.records,
                        idx: &idx[lo_pos..hi],
                    }
                }
            }
            1 | 2 => {
                if first_spanning == lo_pos {
                    if lo_pos == hi {
                        return RecordSet::Empty;
                    }
                    return RecordSet::IndexSlice {
                        records: &self.records,
                        idx: &idx[lo_pos..hi],
                    };
                }

                let prefix_idx = &idx[first_spanning..lo_pos];
                let mut prefix_len = 0usize;
                for &i in prefix_idx {
                    let rec = &self.records[i as usize];
                    if rec.pos <= end && rec.ref_end() >= start {
                        prefix_len += 1;
                    }
                }
                // Records with pos in [start, end]: ref_end >= pos >= start, always overlap.
                let tail = &idx[lo_pos..hi];
                if prefix_len == 0 {
                    if tail.is_empty() {
                        RecordSet::Empty
                    } else {
                        RecordSet::IndexSlice {
                            records: &self.records,
                            idx: tail,
                        }
                    }
                } else {
                    RecordSet::IndexFilteredPrefixAndSlice {
                        records: &self.records,
                        prefix_idx,
                        idx: tail,
                        start,
                        end,
                        prefix_len,
                    }
                }
            }
            _ => RecordSet::Empty,
        }
    }

    pub fn plan_query(
        &self,
        chr: &str,
        start: i64,
        end: i64,
        overlap: u8,
        opts: PlanOptions,
    ) -> (Vec<&VcfRecord>, RegionPlan) {
        let records = self.query(chr, start, end, overlap);
        let plan = plan_region(&records, opts);
        (records, plan)
    }

    pub fn plan_query_set(
        &self,
        chr: &str,
        start: i64,
        end: i64,
        overlap: u8,
        opts: PlanOptions,
    ) -> (RecordSet<'_>, RegionPlan) {
        let records = self.query_set(chr, start, end, overlap);
        let plan = plan_region_set(&records, opts);
        (records, plan)
    }

    fn default_cache_path(path: &Path) -> PathBuf {
        let mut p = path.as_os_str().to_os_string();
        p.push(".cvcf");
        PathBuf::from(p)
    }

    fn tmp_cache_path(path: &Path) -> PathBuf {
        let mut p = path.as_os_str().to_os_string();
        p.push(".tmp");
        PathBuf::from(p)
    }

    fn try_load_default_cache(path: &Path) -> Option<Self> {
        let fp = source_fingerprint(path).ok()?;
        let cache_path = Self::default_cache_path(path);
        Self::read_cache_file(path.to_path_buf(), &cache_path, fp).ok()
    }

    fn write_default_cache(&self) -> io::Result<()> {
        let fp = source_fingerprint(&self.path)?;
        let cache_path = Self::default_cache_path(&self.path);
        let tmp_path = Self::tmp_cache_path(&cache_path);
        self.write_cache_file(&tmp_path, fp)?;
        fs::rename(tmp_path, cache_path)
    }

    fn read_cache_file(
        path: PathBuf,
        cache_path: &Path,
        source_fp: SourceFingerprint,
    ) -> io::Result<Self> {
        let mut r = BufReader::new(File::open(cache_path)?);
        let mut magic = [0u8; 8];
        r.read_exact(&mut magic)?;
        if &magic != CVCF_MAGIC {
            return Err(invalid_data("bad cvcf magic"));
        }
        let version = read_u32(&mut r)?;
        if version != CVCF_VERSION {
            return Err(invalid_data("unsupported cvcf version"));
        }
        let cached_fp = SourceFingerprint {
            len: read_u64(&mut r)?,
            mtime_secs: read_i64(&mut r)?,
            mtime_nanos: read_u32(&mut r)?,
        };
        if cached_fp != source_fp {
            return Err(invalid_data("stale cvcf source fingerprint"));
        }

        let mut store = VcfStore::empty(path);
        let n_seq = read_len(&mut r)?;
        store.seq_names.reserve(n_seq);
        for _ in 0..n_seq {
            let name = read_string(&mut r)?;
            let rid = read_i32(&mut r)?;
            store.seq_names.insert(name, rid);
        }

        let n_samples = read_len(&mut r)?;
        store.sample_names.reserve(n_samples);
        store.sample_idx.reserve(n_samples);
        for i in 0..n_samples {
            let name = read_string(&mut r)?;
            store.sample_idx.insert(name.clone(), i as i32);
            store.sample_names.push(name);
        }
        store.has_gt = read_bool(&mut r)?;
        store.n_sample = read_i32(&mut r)?;

        let n_records = read_len64(&mut r)?;
        store.records.reserve(n_records);
        for _ in 0..n_records {
            let pos = read_i64(&mut r)?;
            let rlen = read_i32(&mut r)?;
            let rid = read_i32(&mut r)?;
            let var_type = read_i32(&mut r)?;

            let n_allele = read_len(&mut r)?;
            let mut alleles: Vec<SmallVec<[u8; 16]>> = Vec::with_capacity(n_allele);
            for _ in 0..n_allele {
                let bytes = read_bytes(&mut r)?;
                alleles.push(SmallVec::from_slice(&bytes));
            }
            let compiled = read_compiled_record(&mut r)?;

            let n_gt_samples = read_len(&mut r)?;
            let mut gt: Vec<SmallVec<[GtAllele; 2]>> = Vec::with_capacity(n_gt_samples);
            for _ in 0..n_gt_samples {
                let n_gt = read_len(&mut r)?;
                let mut sample_gt: SmallVec<[GtAllele; 2]> = SmallVec::new();
                for _ in 0..n_gt {
                    let has_allele = read_bool(&mut r)?;
                    let allele = if has_allele {
                        Some(read_i32(&mut r)?)
                    } else {
                        None
                    };
                    let phased = read_bool(&mut r)?;
                    let raw = read_i32(&mut r)?;
                    sample_gt.push(GtAllele {
                        allele,
                        phased,
                        raw,
                    });
                }
                gt.push(sample_gt);
            }

            let gt_compact = read_compact_gt(&mut r)?;
            let gt_bits = read_gt_bits(&mut r)?;
            store.compile_stats.observe_record(&compiled);
            let (has_gt, is_biallelic_phased_diploid, has_missing_gt) =
                gt_compile_stats_from_stores(n_allele, &gt, gt_compact.as_ref());
            store.compile_stats.observe_gt(
                has_gt,
                is_biallelic_phased_diploid,
                gt_compact.is_some(),
                gt_bits.is_some(),
                has_missing_gt,
            );
            store.has_gt |= has_gt;

            store.records.push(VcfRecord {
                pos,
                rlen,
                rid,
                alleles,
                gt,
                gt_compact,
                gt_bits,
                var_type,
                compiled,
            });
        }
        read_coord_index(&mut r, &mut store)?;
        Ok(store)
    }

    fn write_cache_file(&self, cache_path: &Path, source_fp: SourceFingerprint) -> io::Result<()> {
        let mut w = BufWriter::new(File::create(cache_path)?);
        w.write_all(CVCF_MAGIC)?;
        write_u32(&mut w, CVCF_VERSION)?;
        write_u64(&mut w, source_fp.len)?;
        write_i64(&mut w, source_fp.mtime_secs)?;
        write_u32(&mut w, source_fp.mtime_nanos)?;

        write_len(&mut w, self.seq_names.len())?;
        for (name, rid) in &self.seq_names {
            write_string(&mut w, name)?;
            write_i32(&mut w, *rid)?;
        }

        write_len(&mut w, self.sample_names.len())?;
        for name in &self.sample_names {
            write_string(&mut w, name)?;
        }
        write_bool(&mut w, self.has_gt)?;
        write_i32(&mut w, self.n_sample)?;

        write_len64(&mut w, self.records.len())?;
        for rec in &self.records {
            write_i64(&mut w, rec.pos)?;
            write_i32(&mut w, rec.rlen)?;
            write_i32(&mut w, rec.rid)?;
            write_i32(&mut w, rec.var_type)?;

            write_len(&mut w, rec.alleles.len())?;
            for allele in &rec.alleles {
                write_bytes(&mut w, allele)?;
            }
            write_compiled_record(&mut w, &rec.compiled)?;

            write_len(&mut w, rec.gt.len())?;
            for sample_gt in &rec.gt {
                write_len(&mut w, sample_gt.len())?;
                for gt in sample_gt {
                    write_bool(&mut w, gt.allele.is_some())?;
                    if let Some(allele) = gt.allele {
                        write_i32(&mut w, allele)?;
                    }
                    write_bool(&mut w, gt.phased)?;
                    write_i32(&mut w, gt.raw)?;
                }
            }
            write_compact_gt(&mut w, rec.gt_compact.as_ref())?;
            write_gt_bits(&mut w, rec.gt_bits.as_ref())?;
        }
        write_coord_index(&mut w, self)?;
        w.flush()
    }

    fn rebuild_pmax_end(&mut self) {
        self.pmax_end.clear();
        for (rid, idx) in self.by_rid.iter_mut() {
            idx.sort_by_key(|&i| self.records[i as usize].pos);
            let mut pmax: Vec<i64> = Vec::with_capacity(idx.len());
            let mut m: i64 = i64::MIN;
            for &i in idx.iter() {
                let re = self.records[i as usize].ref_end();
                if re > m {
                    m = re;
                }
                pmax.push(m);
            }
            self.pmax_end.insert(*rid, pmax);
        }
    }

    // -----------------------------------------------------------------------
    // preprocessing
    // -----------------------------------------------------------------------

    fn parse(&mut self) -> Result<(), String> {
        let cpath = CString::new(self.path.to_str().ok_or("non-UTF8 vcf path")?)
            .map_err(|_| "non-NUL vcf path".to_string())?;
        let cmode = CString::new("r").unwrap();

        unsafe {
            let fp = ffi::hts_open(cpath.as_ptr(), cmode.as_ptr());
            if fp.is_null() {
                return Err(format!("hts_open failed for {}", self.path.display()));
            }
            let hdr = ffi::bcf_hdr_read(fp);
            if hdr.is_null() {
                ffi::hts_close(fp);
                return Err(format!("bcf_hdr_read failed for {}", self.path.display()));
            }

            // Samples
            self.n_sample = ffi::shim_bcf_hdr_nsamples(hdr);
            for i in 0..self.n_sample {
                let name_ptr = ffi::shim_bcf_hdr_sample_name(hdr, i);
                let name = if name_ptr.is_null() {
                    format!("sample_{}", i)
                } else {
                    cstr_to_string(name_ptr)
                };
                self.sample_idx.insert(name.clone(), i);
                self.sample_names.push(name);
            }

            // Contig name -> rid map
            let mut nseq: c_int = 0;
            let seqs = ffi::bcf_hdr_seqnames(hdr, &mut nseq);
            if !seqs.is_null() {
                for i in 0..nseq as isize {
                    let p = *seqs.offset(i);
                    if !p.is_null() {
                        let name = cstr_to_string(p);
                        // rid here is the index returned by bcf_hdr_name2id;
                        // bcf_hdr_seqnames indices may differ from rid, so look up.
                        let cname = CString::new(name.as_str()).unwrap();
                        let rid = ffi::shim_bcf_hdr_name2id(hdr, cname.as_ptr());
                        if rid >= 0 {
                            self.seq_names.insert(name, rid);
                        }
                    }
                }
            }

            let rec = ffi::bcf_init();
            if rec.is_null() {
                ffi::bcf_hdr_destroy(hdr);
                ffi::hts_close(fp);
                return Err("bcf_init failed".to_string());
            }

            // GT scratch buffer, grown/reused by htslib across records.
            let mut gt_buf: *mut c_void = std::ptr::null_mut();
            let mut gt_cap: c_int = 0;
            let gt_tag = CString::new("GT").unwrap();
            let unpack_what = ffi::BCF_UN_STR | ffi::BCF_UN_FMT;

            loop {
                let r = ffi::bcf_read(fp, hdr, rec);
                if r == -1 {
                    break; // EOF
                }
                if r < -1 {
                    return Err(format!("bcf_read error code {}", r));
                }
                if ffi::bcf_unpack(rec, unpack_what) < 0 {
                    return Err("bcf_unpack failed".to_string());
                }

                let pos = ffi::shim_bcf_pos(rec);
                let rlen = ffi::shim_bcf_rlen(rec) as i32;
                let rid = ffi::shim_bcf_rid(rec);
                let n_allele = ffi::shim_bcf_n_allele(rec);

                let mut alleles: Vec<SmallVec<[u8; 16]>> = Vec::with_capacity(n_allele as usize);
                for a in 0..n_allele {
                    let p = ffi::shim_bcf_allele(rec, a);
                    let s = if p.is_null() {
                        SmallVec::new()
                    } else {
                        cstr_to_bytes(p)
                    };
                    alleles.push(s);
                }

                // GT
                let ngt = ffi::bcf_get_format_values(
                    hdr,
                    rec,
                    gt_tag.as_ptr(),
                    &mut gt_buf,
                    &mut gt_cap,
                    ffi::BCF_HT_INT,
                );
                let decoded_gt = decode_gt_stores(
                    n_allele as usize,
                    self.n_sample as usize,
                    ngt,
                    gt_buf.cast::<i32>(),
                );
                self.has_gt |= decoded_gt.has_gt;

                let var_type = ffi::bcf_get_variant_types(rec);
                let compiled = CompiledRecord::from_alleles(rlen, &alleles);
                self.compile_stats.observe_record(&compiled);
                self.compile_stats.observe_gt(
                    decoded_gt.has_gt,
                    decoded_gt.is_biallelic_phased_diploid,
                    decoded_gt.gt_compact.is_some(),
                    decoded_gt.gt_bits.is_some(),
                    decoded_gt.has_missing_gt,
                );

                let rid_bucket = self.by_rid.entry(rid).or_default();
                rid_bucket.push(self.records.len() as u32);
                self.records.push(VcfRecord {
                    pos,
                    rlen,
                    rid,
                    alleles,
                    gt: decoded_gt.gt,
                    gt_compact: decoded_gt.gt_compact,
                    gt_bits: decoded_gt.gt_bits,
                    var_type,
                    compiled,
                });
            }

            if !gt_buf.is_null() {
                ffi::free(gt_buf);
            }
            ffi::bcf_destroy(rec);
            ffi::bcf_hdr_destroy(hdr);
            ffi::hts_close(fp);
        }

        // Sort each rid bucket by pos (stable on file order) and build pmax_end.
        self.rebuild_pmax_end();

        Ok(())
    }
}

// ---------------------------------------------------------------------------
// compact owned-file cache helpers
// ---------------------------------------------------------------------------

fn source_fingerprint(path: &Path) -> io::Result<SourceFingerprint> {
    let meta = fs::metadata(path)?;
    let modified = meta.modified()?;
    let duration = modified
        .duration_since(UNIX_EPOCH)
        .map_err(|_| invalid_data("source mtime before unix epoch"))?;
    Ok(SourceFingerprint {
        len: meta.len(),
        mtime_secs: duration.as_secs() as i64,
        mtime_nanos: duration.subsec_nanos(),
    })
}

struct DecodedGtStores {
    gt: Vec<SmallVec<[GtAllele; 2]>>,
    gt_compact: Option<CompactGt>,
    gt_bits: Option<BiallelicPhasedGtBits>,
    has_gt: bool,
    is_biallelic_phased_diploid: bool,
    has_missing_gt: bool,
}

fn decode_gt_stores(
    n_allele: usize,
    n_samples: usize,
    ngt: c_int,
    base: *const i32,
) -> DecodedGtStores {
    if ngt <= 0 || n_samples == 0 || base.is_null() {
        return DecodedGtStores {
            gt: Vec::new(),
            gt_compact: None,
            gt_bits: None,
            has_gt: false,
            is_biallelic_phased_diploid: false,
            has_missing_gt: false,
        };
    }

    // bcf_get_format_values returns nsmpl*max_ploidy int32 values.
    let ploidy = (ngt as usize) / n_samples;
    let mut compact_offsets = Vec::with_capacity(n_samples + 1);
    let mut compact_codes = Vec::with_capacity(ngt as usize);
    compact_offsets.push(0);
    let mut compact_ok = true;
    let mut has_missing_gt = false;
    let mut is_biallelic_phased_diploid = n_allele == 2;
    let mut gt_bits = (n_allele == 2).then(|| BiallelicPhasedGtBits {
        n_samples,
        hap1_alt: vec![0; n_samples.div_ceil(64)],
        hap2_alt: vec![0; n_samples.div_ceil(64)],
        missing: vec![0; n_samples.div_ceil(64)],
    });

    for sample_idx in 0..n_samples {
        let mut sample_len = 0usize;
        let mut hap_raw = [0i32; 2];
        for j in 0..ploidy {
            let raw = unsafe { *base.add(sample_idx * ploidy + j) };
            if raw == ffi::BCF_INT32_VECTOR_END {
                break;
            }
            if sample_len < 2 {
                hap_raw[sample_len] = raw;
            }
            sample_len += 1;

            if ffi::gt_is_missing(raw) {
                has_missing_gt = true;
            }
            if compact_ok {
                match compact_gt_code(raw) {
                    Some(code) => compact_codes.push(code),
                    None => {
                        compact_ok = false;
                        compact_offsets.clear();
                        compact_codes.clear();
                    }
                }
            }
        }

        if compact_ok {
            if let Ok(offset) = u32::try_from(compact_codes.len()) {
                compact_offsets.push(offset);
            } else {
                compact_ok = false;
                compact_offsets.clear();
                compact_codes.clear();
            }
        }

        if is_biallelic_phased_diploid && sample_len != 2 {
            is_biallelic_phased_diploid = false;
        }
        if let Some(bits) = gt_bits.as_mut() {
            if sample_len != 2 {
                gt_bits = None;
                continue;
            }
            let word = sample_idx / 64;
            let bit = 1u64 << (sample_idx & 63);
            let a0_missing = ffi::gt_is_missing(hap_raw[0]);
            let a1_missing = ffi::gt_is_missing(hap_raw[1]);
            if a0_missing || a1_missing {
                bits.missing[word] |= bit;
                is_biallelic_phased_diploid = false;
                continue;
            }
            let h0 = ffi::gt_allele(hap_raw[0]);
            let h1 = ffi::gt_allele(hap_raw[1]);
            if !ffi::gt_is_phased(hap_raw[0])
                || !ffi::gt_is_phased(hap_raw[1])
                || !(0..=1).contains(&h0)
                || !(0..=1).contains(&h1)
            {
                gt_bits = None;
                is_biallelic_phased_diploid = false;
                continue;
            }
            if h0 == 1 {
                bits.hap1_alt[word] |= bit;
            }
            if h1 == 1 {
                bits.hap2_alt[word] |= bit;
            }
        }
    }

    let gt_compact = compact_ok.then_some(CompactGt {
        sample_offsets: compact_offsets,
        codes: compact_codes,
    });
    let gt = if gt_compact.is_some() {
        Vec::new()
    } else {
        decode_raw_gt_matrix(n_samples, ngt, base)
    };

    DecodedGtStores {
        gt,
        gt_compact,
        gt_bits,
        has_gt: true,
        is_biallelic_phased_diploid,
        has_missing_gt,
    }
}

#[inline]
fn compact_gt_code(raw: i32) -> Option<u16> {
    let mut code = if ffi::gt_is_missing(raw) {
        COMPACT_GT_MISSING
    } else {
        let idx = ffi::gt_allele(raw);
        if idx < 0 || idx > COMPACT_GT_ALLELE_MASK as i32 {
            return None;
        }
        idx as u16
    };
    if ffi::gt_is_phased(raw) {
        code |= COMPACT_GT_PHASED;
    }
    Some(code)
}

fn decode_raw_gt_matrix(
    n_samples: usize,
    ngt: c_int,
    base: *const i32,
) -> Vec<SmallVec<[GtAllele; 2]>> {
    let ploidy = (ngt as usize) / n_samples;
    let mut gt = Vec::with_capacity(n_samples);
    for sample_idx in 0..n_samples {
        let mut sample_gt: SmallVec<[GtAllele; 2]> = SmallVec::new();
        for j in 0..ploidy {
            let raw = unsafe { *base.add(sample_idx * ploidy + j) };
            if raw == ffi::BCF_INT32_VECTOR_END {
                break;
            }
            let allele = if ffi::gt_is_missing(raw) {
                None
            } else {
                Some(ffi::gt_allele(raw))
            };
            sample_gt.push(GtAllele {
                allele,
                phased: ffi::gt_is_phased(raw),
                raw,
            });
        }
        gt.push(sample_gt);
    }
    gt
}

fn gt_compile_stats_from_stores(
    n_allele: usize,
    gt: &[SmallVec<[GtAllele; 2]>],
    compact: Option<&CompactGt>,
) -> (bool, bool, bool) {
    if !gt.is_empty() {
        return gt_compile_stats(n_allele, gt);
    }
    let Some(compact) = compact else {
        return (false, false, false);
    };
    let mut has_missing = false;
    let mut is_biallelic_phased_diploid = n_allele == 2;
    for sample_idx in 0..compact.n_samples() {
        let Some(sample_gt) = compact.sample(sample_idx) else {
            is_biallelic_phased_diploid = false;
            continue;
        };
        if sample_gt.len() != 2 {
            is_biallelic_phased_diploid = false;
        }
        for &code in sample_gt {
            match CompactGt::allele(code) {
                Some(allele) => {
                    if !CompactGt::phased(code) || (n_allele == 2 && !(0..=1).contains(&allele)) {
                        is_biallelic_phased_diploid = false;
                    }
                }
                None => {
                    has_missing = true;
                    is_biallelic_phased_diploid = false;
                }
            }
        }
    }
    (true, is_biallelic_phased_diploid, has_missing)
}

fn gt_compile_stats(n_allele: usize, gt: &[SmallVec<[GtAllele; 2]>]) -> (bool, bool, bool) {
    if gt.is_empty() {
        return (false, false, false);
    }
    let mut has_missing = false;
    let mut is_biallelic_phased_diploid = n_allele == 2;
    for sample_gt in gt {
        if sample_gt.len() != 2 {
            is_biallelic_phased_diploid = false;
        }
        for allele in sample_gt {
            if allele.allele.is_none() {
                has_missing = true;
                is_biallelic_phased_diploid = false;
            }
            if !allele.phased {
                is_biallelic_phased_diploid = false;
            }
        }
    }
    (true, is_biallelic_phased_diploid, has_missing)
}

fn invalid_data(msg: &'static str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg)
}

fn checked_u32_len(n: usize) -> io::Result<u32> {
    u32::try_from(n).map_err(|_| invalid_data("length exceeds u32"))
}

fn checked_usize_len(n: u64) -> io::Result<usize> {
    usize::try_from(n).map_err(|_| invalid_data("length exceeds usize"))
}

fn read_len<R: Read>(r: &mut R) -> io::Result<usize> {
    Ok(read_u32(r)? as usize)
}

fn read_len64<R: Read>(r: &mut R) -> io::Result<usize> {
    checked_usize_len(read_u64(r)?)
}

fn write_len<W: Write>(w: &mut W, n: usize) -> io::Result<()> {
    write_u32(w, checked_u32_len(n)?)
}

fn write_len64<W: Write>(w: &mut W, n: usize) -> io::Result<()> {
    write_u64(w, n as u64)
}

fn read_bool<R: Read>(r: &mut R) -> io::Result<bool> {
    match read_u8(r)? {
        0 => Ok(false),
        1 => Ok(true),
        _ => Err(invalid_data("invalid bool")),
    }
}

fn write_bool<W: Write>(w: &mut W, v: bool) -> io::Result<()> {
    write_u8(w, u8::from(v))
}

fn read_string<R: Read>(r: &mut R) -> io::Result<String> {
    let bytes = read_bytes(r)?;
    String::from_utf8(bytes).map_err(|_| invalid_data("invalid utf8 string"))
}

fn write_string<W: Write>(w: &mut W, s: &str) -> io::Result<()> {
    write_bytes(w, s.as_bytes())
}

fn read_bytes<R: Read>(r: &mut R) -> io::Result<Vec<u8>> {
    let len = read_len(r)?;
    if len > (1usize << 30) {
        return Err(invalid_data("byte field too large"));
    }
    let mut bytes = vec![0u8; len];
    r.read_exact(&mut bytes)?;
    Ok(bytes)
}

fn write_bytes<W: Write>(w: &mut W, bytes: &[u8]) -> io::Result<()> {
    write_len(w, bytes.len())?;
    w.write_all(bytes)
}

fn read_compiled_record<R: Read>(r: &mut R) -> io::Result<CompiledRecord> {
    let kind = read_record_kind(r)?;
    let flags = RecordFlags::from_bits(read_u16(r)?);
    let n_ops = read_len(r)?;
    if n_ops > (1usize << 20) {
        return Err(invalid_data("compiled op table too large"));
    }
    let mut ops: SmallVec<[AlleleOp; 2]> = SmallVec::with_capacity(n_ops);
    for _ in 0..n_ops {
        ops.push(AlleleOp {
            kind: read_allele_op_kind(r)?,
            ref_len: read_u32(r)?,
            alt_len: read_u32(r)?,
            trim_beg: read_u16(r)?,
            len_diff: read_i32(r)?,
            case_flags: read_u8(r)?,
        });
    }
    Ok(CompiledRecord { kind, flags, ops })
}

fn write_compiled_record<W: Write>(w: &mut W, compiled: &CompiledRecord) -> io::Result<()> {
    write_u8(w, compiled.kind as u8)?;
    write_u16(w, compiled.flags.bits())?;
    write_len(w, compiled.ops.len())?;
    for op in &compiled.ops {
        write_u8(w, op.kind as u8)?;
        write_u32(w, op.ref_len)?;
        write_u32(w, op.alt_len)?;
        write_u16(w, op.trim_beg)?;
        write_i32(w, op.len_diff)?;
        write_u8(w, op.case_flags)?;
    }
    Ok(())
}

fn read_record_kind<R: Read>(r: &mut R) -> io::Result<RecordKind> {
    match read_u8(r)? {
        0 => Ok(RecordKind::RefOnly),
        1 => Ok(RecordKind::Snp1),
        2 => Ok(RecordKind::SameLen),
        3 => Ok(RecordKind::NormInsertion),
        4 => Ok(RecordKind::NormDeletion),
        5 => Ok(RecordKind::SimpleIndel),
        6 => Ok(RecordKind::SymbolicDel),
        7 => Ok(RecordKind::GvcfBlock),
        8 => Ok(RecordKind::Complex),
        _ => Err(invalid_data("invalid record kind")),
    }
}

fn read_allele_op_kind<R: Read>(r: &mut R) -> io::Result<AlleleOpKind> {
    match read_u8(r)? {
        0 => Ok(AlleleOpKind::Ref),
        1 => Ok(AlleleOpKind::SameLen),
        2 => Ok(AlleleOpKind::Insert),
        3 => Ok(AlleleOpKind::Delete),
        4 => Ok(AlleleOpKind::Replace),
        5 => Ok(AlleleOpKind::SymbolicDel),
        6 => Ok(AlleleOpKind::GvcfRefBlock),
        7 => Ok(AlleleOpKind::Missing),
        8 => Ok(AlleleOpKind::Unsupported),
        _ => Err(invalid_data("invalid allele op kind")),
    }
}

fn read_compact_gt<R: Read>(r: &mut R) -> io::Result<Option<CompactGt>> {
    if !read_bool(r)? {
        return Ok(None);
    }
    let sample_offsets = read_u32_vec(r)?;
    let codes = read_u16_vec(r)?;
    if sample_offsets.is_empty() || sample_offsets[0] != 0 {
        return Err(invalid_data("invalid compact GT offsets"));
    }
    let mut prev = 0u32;
    for &offset in &sample_offsets {
        if offset < prev || offset as usize > codes.len() {
            return Err(invalid_data("invalid compact GT offset order"));
        }
        prev = offset;
    }
    if sample_offsets.last().copied().unwrap_or(0) as usize != codes.len() {
        return Err(invalid_data("compact GT offsets do not cover codes"));
    }
    Ok(Some(CompactGt {
        sample_offsets,
        codes,
    }))
}

fn write_compact_gt<W: Write>(w: &mut W, compact: Option<&CompactGt>) -> io::Result<()> {
    let Some(compact) = compact else {
        return write_bool(w, false);
    };
    write_bool(w, true)?;
    write_u32_slice(w, &compact.sample_offsets)?;
    write_u16_slice(w, &compact.codes)
}

fn read_gt_bits<R: Read>(r: &mut R) -> io::Result<Option<BiallelicPhasedGtBits>> {
    if !read_bool(r)? {
        return Ok(None);
    }
    let n_samples = read_len(r)?;
    let hap1_alt = read_u64_vec(r)?;
    let hap2_alt = read_u64_vec(r)?;
    let missing = read_u64_vec(r)?;
    let expected_words = n_samples.div_ceil(64);
    if hap1_alt.len() != expected_words
        || hap2_alt.len() != expected_words
        || missing.len() != expected_words
    {
        return Err(invalid_data("invalid gt bitset length"));
    }
    Ok(Some(BiallelicPhasedGtBits {
        n_samples,
        hap1_alt,
        hap2_alt,
        missing,
    }))
}

fn write_gt_bits<W: Write>(w: &mut W, bits: Option<&BiallelicPhasedGtBits>) -> io::Result<()> {
    let Some(bits) = bits else {
        return write_bool(w, false);
    };
    write_bool(w, true)?;
    write_len(w, bits.n_samples)?;
    write_u64_slice(w, &bits.hap1_alt)?;
    write_u64_slice(w, &bits.hap2_alt)?;
    write_u64_slice(w, &bits.missing)
}

fn read_coord_index<R: Read>(r: &mut R, store: &mut VcfStore) -> io::Result<()> {
    store.by_rid.clear();
    store.pmax_end.clear();
    let n_buckets = read_len(r)?;
    store.by_rid.reserve(n_buckets);
    store.pmax_end.reserve(n_buckets);
    for _ in 0..n_buckets {
        let rid = read_i32(r)?;
        let idx = read_u32_vec(r)?;
        let pmax = read_i64_vec(r)?;
        if idx.len() != pmax.len() {
            return Err(invalid_data("coord index length mismatch"));
        }
        let mut prev_pos = i64::MIN;
        let mut prev_pmax = i64::MIN;
        for (&record_idx, &pmax_value) in idx.iter().zip(&pmax) {
            let rec = store
                .records
                .get(record_idx as usize)
                .ok_or_else(|| invalid_data("coord index record out of range"))?;
            if rec.rid != rid {
                return Err(invalid_data("coord index rid mismatch"));
            }
            if rec.pos < prev_pos || pmax_value < prev_pmax || pmax_value < rec.ref_end() {
                return Err(invalid_data("invalid coord index ordering"));
            }
            prev_pos = rec.pos;
            prev_pmax = pmax_value;
        }
        store.by_rid.insert(rid, idx);
        store.pmax_end.insert(rid, pmax);
    }
    Ok(())
}

fn write_coord_index<W: Write>(w: &mut W, store: &VcfStore) -> io::Result<()> {
    write_len(w, store.by_rid.len())?;
    for (&rid, idx) in &store.by_rid {
        write_i32(w, rid)?;
        write_u32_slice(w, idx)?;
        let pmax = store
            .pmax_end
            .get(&rid)
            .ok_or_else(|| invalid_data("missing pmax coord index"))?;
        write_i64_slice(w, pmax)?;
    }
    Ok(())
}

fn read_u64_vec<R: Read>(r: &mut R) -> io::Result<Vec<u64>> {
    let len = read_len(r)?;
    if len > (1usize << 28) {
        return Err(invalid_data("u64 vector too large"));
    }
    #[cfg(target_endian = "little")]
    {
        let byte_len = len
            .checked_mul(std::mem::size_of::<u64>())
            .ok_or_else(|| invalid_data("u64 vector byte length overflow"))?;
        let mut out = vec![0u64; len];
        // The cache format is little-endian and this branch only compiles on
        // little-endian targets; the initialized numeric buffer can be viewed as
        // bytes for one contiguous read.
        let bytes =
            unsafe { std::slice::from_raw_parts_mut(out.as_mut_ptr().cast::<u8>(), byte_len) };
        r.read_exact(bytes)?;
        Ok(out)
    }
    #[cfg(not(target_endian = "little"))]
    {
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            out.push(read_u64(r)?);
        }
        Ok(out)
    }
}

fn write_u64_slice<W: Write>(w: &mut W, xs: &[u64]) -> io::Result<()> {
    write_len(w, xs.len())?;
    #[cfg(target_endian = "little")]
    {
        let byte_len = xs
            .len()
            .checked_mul(std::mem::size_of::<u64>())
            .ok_or_else(|| invalid_data("u64 slice byte length overflow"))?;
        let bytes = unsafe { std::slice::from_raw_parts(xs.as_ptr().cast::<u8>(), byte_len) };
        w.write_all(bytes)
    }
    #[cfg(not(target_endian = "little"))]
    {
        for &x in xs {
            write_u64(w, x)?;
        }
        Ok(())
    }
}

fn read_u32_vec<R: Read>(r: &mut R) -> io::Result<Vec<u32>> {
    let len = read_len(r)?;
    if len > (1usize << 30) {
        return Err(invalid_data("u32 vector too large"));
    }
    #[cfg(target_endian = "little")]
    {
        let byte_len = len
            .checked_mul(std::mem::size_of::<u32>())
            .ok_or_else(|| invalid_data("u32 vector byte length overflow"))?;
        let mut out = vec![0u32; len];
        let bytes =
            unsafe { std::slice::from_raw_parts_mut(out.as_mut_ptr().cast::<u8>(), byte_len) };
        r.read_exact(bytes)?;
        Ok(out)
    }
    #[cfg(not(target_endian = "little"))]
    {
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            out.push(read_u32(r)?);
        }
        Ok(out)
    }
}

fn write_u32_slice<W: Write>(w: &mut W, xs: &[u32]) -> io::Result<()> {
    write_len(w, xs.len())?;
    #[cfg(target_endian = "little")]
    {
        let byte_len = xs
            .len()
            .checked_mul(std::mem::size_of::<u32>())
            .ok_or_else(|| invalid_data("u32 slice byte length overflow"))?;
        let bytes = unsafe { std::slice::from_raw_parts(xs.as_ptr().cast::<u8>(), byte_len) };
        w.write_all(bytes)
    }
    #[cfg(not(target_endian = "little"))]
    {
        for &x in xs {
            write_u32(w, x)?;
        }
        Ok(())
    }
}

fn read_u16_vec<R: Read>(r: &mut R) -> io::Result<Vec<u16>> {
    let len = read_len(r)?;
    if len > (1usize << 30) {
        return Err(invalid_data("u16 vector too large"));
    }
    #[cfg(target_endian = "little")]
    {
        let byte_len = len
            .checked_mul(std::mem::size_of::<u16>())
            .ok_or_else(|| invalid_data("u16 vector byte length overflow"))?;
        let mut out = vec![0u16; len];
        let bytes =
            unsafe { std::slice::from_raw_parts_mut(out.as_mut_ptr().cast::<u8>(), byte_len) };
        r.read_exact(bytes)?;
        Ok(out)
    }
    #[cfg(not(target_endian = "little"))]
    {
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            out.push(read_u16(r)?);
        }
        Ok(out)
    }
}

fn write_u16_slice<W: Write>(w: &mut W, xs: &[u16]) -> io::Result<()> {
    write_len(w, xs.len())?;
    #[cfg(target_endian = "little")]
    {
        let byte_len = xs
            .len()
            .checked_mul(std::mem::size_of::<u16>())
            .ok_or_else(|| invalid_data("u16 slice byte length overflow"))?;
        let bytes = unsafe { std::slice::from_raw_parts(xs.as_ptr().cast::<u8>(), byte_len) };
        w.write_all(bytes)
    }
    #[cfg(not(target_endian = "little"))]
    {
        for &x in xs {
            write_u16(w, x)?;
        }
        Ok(())
    }
}

fn read_i64_vec<R: Read>(r: &mut R) -> io::Result<Vec<i64>> {
    let len = read_len(r)?;
    if len > (1usize << 30) {
        return Err(invalid_data("i64 vector too large"));
    }
    #[cfg(target_endian = "little")]
    {
        let byte_len = len
            .checked_mul(std::mem::size_of::<i64>())
            .ok_or_else(|| invalid_data("i64 vector byte length overflow"))?;
        let mut out = vec![0i64; len];
        let bytes =
            unsafe { std::slice::from_raw_parts_mut(out.as_mut_ptr().cast::<u8>(), byte_len) };
        r.read_exact(bytes)?;
        Ok(out)
    }
    #[cfg(not(target_endian = "little"))]
    {
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            out.push(read_i64(r)?);
        }
        Ok(out)
    }
}

fn write_i64_slice<W: Write>(w: &mut W, xs: &[i64]) -> io::Result<()> {
    write_len(w, xs.len())?;
    #[cfg(target_endian = "little")]
    {
        let byte_len = xs
            .len()
            .checked_mul(std::mem::size_of::<i64>())
            .ok_or_else(|| invalid_data("i64 slice byte length overflow"))?;
        let bytes = unsafe { std::slice::from_raw_parts(xs.as_ptr().cast::<u8>(), byte_len) };
        w.write_all(bytes)
    }
    #[cfg(not(target_endian = "little"))]
    {
        for &x in xs {
            write_i64(w, x)?;
        }
        Ok(())
    }
}

fn read_u8<R: Read>(r: &mut R) -> io::Result<u8> {
    let mut b = [0u8; 1];
    r.read_exact(&mut b)?;
    Ok(b[0])
}

fn write_u8<W: Write>(w: &mut W, v: u8) -> io::Result<()> {
    w.write_all(&[v])
}

fn read_u16<R: Read>(r: &mut R) -> io::Result<u16> {
    let mut b = [0u8; 2];
    r.read_exact(&mut b)?;
    Ok(u16::from_le_bytes(b))
}

fn write_u16<W: Write>(w: &mut W, v: u16) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

fn read_u32<R: Read>(r: &mut R) -> io::Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

fn write_u32<W: Write>(w: &mut W, v: u32) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

fn read_i32<R: Read>(r: &mut R) -> io::Result<i32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(i32::from_le_bytes(b))
}

fn write_i32<W: Write>(w: &mut W, v: i32) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

fn read_u64<R: Read>(r: &mut R) -> io::Result<u64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(u64::from_le_bytes(b))
}

fn write_u64<W: Write>(w: &mut W, v: u64) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

fn read_i64<R: Read>(r: &mut R) -> io::Result<i64> {
    let mut b = [0u8; 8];
    r.read_exact(&mut b)?;
    Ok(i64::from_le_bytes(b))
}

fn write_i64<W: Write>(w: &mut W, v: i64) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

// ---------------------------------------------------------------------------
// small C-string helpers
// ---------------------------------------------------------------------------

unsafe fn cstr_to_string(p: *const std::os::raw::c_char) -> String {
    let s = std::ffi::CStr::from_ptr(p);
    s.to_string_lossy().into_owned()
}

unsafe fn cstr_to_bytes(p: *const std::os::raw::c_char) -> SmallVec<[u8; 16]> {
    let s = std::ffi::CStr::from_ptr(p);
    SmallVec::from_slice(s.to_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a tiny VCF (with GT) on disk and return its path.
    /// `name` makes the temp dir unique per test so parallel `cargo test`
    /// doesn't race on a shared directory.
    fn write_vcf(name: &str, body: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("consensus_rs_vcf_{}", name));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let vcf = dir.join("test.vcf");
        let header = "##fileformat=VCFv4.3\n\
            ##contig=<ID=chr1,length=1000>\n\
            ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";
        std::fs::write(&vcf, format!("{}{}", header, body)).unwrap();
        vcf
    }

    #[test]
    fn parses_records_samples_and_gt() {
        // 1-based POS in VCF; htslib gives 0-based pos.
        // chr1:10 SNP, chr1:20 del (REF=AC, ALT=A -> del of 1), chr1:30 ins
        let vcf = write_vcf(
            "parses",
            "chr1\t10\t.\tG\tA\t.\t.\t.\tGT\t0|1\t1/1\n\
             chr1\t20\t.\tAC\tA\t.\t.\t.\tGT\t0/0\t./.\n\
             chr1\t30\t.\tA\tAG\t.\t.\t.\tGT\t1|1\t0|1\n",
        );
        let store = VcfStore::load(&vcf).unwrap();
        assert_eq!(store.n_sample(), 2);
        assert_eq!(store.sample_names(), &["S1", "S2"]);
        assert_eq!(store.n_records(), 3);
        assert!(store.has_gt());
        assert_eq!(store.compile_stats().compact_gt_records, 3);

        let r0 = &store.records()[0];
        assert_eq!(r0.pos, 9); // 1-based 10 -> 0-based 9
        assert_eq!(r0.rlen, 1);
        assert_eq!(r0.alleles.len(), 2);
        assert_eq!(&r0.alleles[0][..], b"G");
        assert_eq!(&r0.alleles[1][..], b"A");
        assert!(
            r0.gt.is_empty(),
            "compact GT records do not keep the raw matrix on the hot path"
        );
        // S1 = 0|1 (phased): allele0=REF(0) phased, allele1=ALT(1) phased
        let compact = r0.gt_compact.as_ref().expect("compact GT store");
        let s0 = compact.sample(0).unwrap();
        assert_eq!(s0.len(), 2);
        assert_eq!(CompactGt::allele(s0[0]), Some(0));
        assert!(CompactGt::phased(s0[0]));
        assert_eq!(CompactGt::allele(s0[1]), Some(1));
        assert!(CompactGt::phased(s0[1]));
        // S2 = 1/1 (unphased): both ALT, unphased
        let s1 = compact.sample(1).unwrap();
        assert_eq!(CompactGt::allele(s1[0]), Some(1));
        assert!(!CompactGt::phased(s1[0]));

        // missing GT ./.
        let r1 = &store.records()[1];
        let compact = r1.gt_compact.as_ref().expect("compact GT store");
        assert_eq!(CompactGt::allele(compact.sample(1).unwrap()[0]), None);
    }

    #[test]
    fn writes_and_reads_owned_cvcf_cache() {
        let vcf = write_vcf(
            "cache_roundtrip",
            "chr1\t10\t.\tG\tA\t.\t.\t.\tGT\t0|1\t1|1\n\
             chr1\t20\t.\tAC\tA\t.\t.\t.\tGT\t0/0\t./.\n",
        );
        let cache_path = VcfStore::default_cache_path(&vcf);
        assert!(!cache_path.exists());

        let parsed = VcfStore::load(&vcf).unwrap();
        assert!(cache_path.exists());

        let fp = source_fingerprint(&vcf).unwrap();
        let cached = VcfStore::read_cache_file(vcf.clone(), &cache_path, fp).unwrap();

        assert_eq!(cached.n_records(), parsed.n_records());
        assert_eq!(cached.sample_names(), parsed.sample_names());
        assert_eq!(
            cached.compile_stats().records_total,
            parsed.compile_stats().records_total
        );
        assert_eq!(
            cached.compile_stats().compact_gt_records,
            parsed.compile_stats().compact_gt_records
        );
        let q = cached.query("chr1", 0, 25, 1);
        assert_eq!(q.len(), 2);
        assert_eq!(&q[0].alleles[1][..], b"A");
        assert!(q[1].gt.is_empty());
        assert_eq!(
            CompactGt::allele(q[1].gt_compact.as_ref().unwrap().sample(1).unwrap()[0]),
            None
        );
        let spanning = cached.query("chr1", 20, 20, 1);
        assert_eq!(spanning.len(), 1);
        assert_eq!(spanning[0].pos, 19);
        assert!(cached.query("chr1", 20, 20, 0).is_empty());
        assert_eq!(cached.records()[0].compiled, parsed.records()[0].compiled);
        assert_eq!(cached.records()[1].compiled, parsed.records()[1].compiled);
        let compact = cached.records()[0]
            .gt_compact
            .as_ref()
            .expect("cache should persist compact GT");
        let s1 = compact.sample(1).unwrap();
        assert_eq!(CompactGt::allele(s1[0]), Some(1));
        assert_eq!(CompactGt::allele(s1[1]), Some(1));
        let bits = q[0]
            .gt_bits
            .as_ref()
            .expect("cache should persist precompiled GT bitset");
        assert_eq!(bits.allele_for_hap(0, 1), Some(Some(0)));
        assert_eq!(bits.allele_for_hap(0, 2), Some(Some(1)));
    }

    #[test]
    fn builds_biallelic_phased_gt_bitset() {
        let vcf = write_vcf("gt_bitset", "chr1\t10\t.\tG\tA\t.\t.\t.\tGT\t0|1\t1|0\n");
        let store = VcfStore::load(&vcf).unwrap();
        let rec = &store.records()[0];
        let bits = rec.gt_bits.as_ref().expect("biallelic phased bitset");

        assert_eq!(bits.allele_for_hap(0, 1), Some(Some(0)));
        assert_eq!(bits.allele_for_hap(0, 2), Some(Some(1)));
        assert_eq!(bits.allele_for_hap(1, 1), Some(Some(1)));
        assert_eq!(bits.allele_for_hap(1, 2), Some(Some(0)));
        assert_eq!(store.compile_stats().biallelic_gt_bitset_records, 1);
    }

    #[test]
    fn decode_gt_keeps_raw_matrix_only_when_compact_overflows() {
        let compact_raw = [((1 + 1) << 1) | 1];
        let compact = decode_gt_stores(2, 1, compact_raw.len() as c_int, compact_raw.as_ptr());
        assert!(compact.has_gt);
        assert!(compact.gt.is_empty());
        assert!(compact.gt_compact.is_some());

        let huge_allele = COMPACT_GT_ALLELE_MASK as i32 + 1;
        let raw = [(huge_allele + 1) << 1];
        let fallback = decode_gt_stores(
            huge_allele as usize + 1,
            1,
            raw.len() as c_int,
            raw.as_ptr(),
        );
        assert!(fallback.has_gt);
        assert!(fallback.gt_compact.is_none());
        assert_eq!(fallback.gt.len(), 1);
        assert_eq!(fallback.gt[0][0].allele, Some(huge_allele));
    }

    #[test]
    fn region_query_does_not_miss_spanning_deletion() {
        // A deletion starting BEFORE the region but spanning into it.
        // chr1:5 REF=GGGGG ALT=G (deletes 4 bases spanning 1-based 6..9)
        // plus a SNP inside the region at chr1:8
        let vcf = write_vcf(
            "spanning_del",
            "chr1\t5\t.\tGGGGG\tG\t.\t.\t.\tGT\t1|1\t0/1\n\
             chr1\t8\t.\tC\tT\t.\t.\t.\tGT\t0/1\t1/1\n",
        );
        let store = VcfStore::load(&vcf).unwrap();

        // region 0-based [6, 12] (1-based 7..13). The deletion pos=4 (1-based 5),
        // ref_end = 4 + 5 - 1 = 8 >= 6 -> must be included for overlap=1/2.
        let q1 = store.query("chr1", 6, 12, 1);
        assert_eq!(q1.len(), 2, "overlap=1 must include the spanning deletion");
        assert_eq!(q1[0].pos, 4);
        assert_eq!(q1[1].pos, 7);

        let q2 = store.query("chr1", 6, 12, 2);
        assert_eq!(q2.len(), 2, "overlap=2 must include the spanning deletion");

        // overlap=0 (POS in region) excludes the deletion (pos=4 < 6)
        let q0 = store.query("chr1", 6, 12, 0);
        assert_eq!(q0.len(), 1);
        assert_eq!(q0[0].pos, 7);
    }

    #[test]
    fn query_set_borrows_hot_path_and_filtered_spanning_prefix() {
        let vcf = write_vcf(
            "query_set_shape",
            "chr1\t5\t.\tGGGGG\tG\t.\t.\t.\tGT\t1|1\t0/1\n\
             chr1\t8\t.\tC\tT\t.\t.\t.\tGT\t0/1\t1/1\n\
             chr1\t12\t.\tA\tG\t.\t.\t.\tGT\t0/1\t1/1\n",
        );
        let store = VcfStore::load(&vcf).unwrap();

        let hot = store.query_set("chr1", 9, 20, 0);
        match &hot {
            RecordSet::IndexSlice { idx, .. } => assert_eq!(idx.len(), 1),
            other => panic!("expected borrowed hot-path slice, got {:?}", other.len()),
        }
        assert_eq!(hot.iter().next().unwrap().pos, 11);

        let spanning = store.query_set("chr1", 6, 12, 1);
        match &spanning {
            RecordSet::IndexFilteredPrefixAndSlice {
                prefix_len, idx, ..
            } => {
                assert_eq!(*prefix_len, 1);
                assert_eq!(idx.len(), 2);
            }
            other => panic!(
                "expected filtered spanning prefix + borrowed tail, got {:?}",
                other.len()
            ),
        }
        let positions: Vec<i64> = spanning.iter().map(|r| r.pos).collect();
        assert_eq!(positions, vec![4, 7, 11]);
    }

    #[test]
    fn region_query_overlap0_pos_in_region() {
        let vcf = write_vcf(
            "overlap0",
            "chr1\t10\t.\tG\tA\t.\t.\t.\tGT\t0/1\t0/1\n\
             chr1\t20\t.\tG\tA\t.\t.\t.\tGT\t0/1\t0/1\n\
             chr1\t30\t.\tG\tA\t.\t.\t.\tGT\t0/1\t0/1\n",
        );
        let store = VcfStore::load(&vcf).unwrap();
        // 0-based [9, 28] -> 1-based 10..29 -> records at pos 9 and 19
        let q = store.query("chr1", 9, 28, 0);
        assert_eq!(q.len(), 2);
        assert_eq!(q[0].pos, 9);
        assert_eq!(q[1].pos, 19);
    }

    /// Parity with `bcftools view -r` record count (M1 acceptance). Ignored by
    /// default so `cargo test` stays hermetic; run with `-- --ignored` and
    /// bgzip/tabix/bcftools on PATH.
    #[test]
    #[ignore]
    fn bcftools_view_region_count_parity() {
        use std::process::Command;
        let vcf = write_vcf(
            "bcftools_parity",
            // pos (1-based): 5 del(GGGGG>G), 12 snp, 18 mnp(AC>TA), 25 ins, 30 snp
            "chr1\t5\t.\tGGGGG\tG\t.\t.\t.\tGT\t0|1\t1/1\n\
             chr1\t12\t.\tG\tA\t.\t.\t.\tGT\t0/1\t0/1\n\
             chr1\t18\t.\tAC\tTA\t.\t.\t.\tGT\t1|0\t./.\n\
             chr1\t25\t.\tA\tAGT\t.\t.\t.\tGT\t0/1\t1/1\n\
             chr1\t30\t.\tC\tG\t.\t.\t.\tGT\t0/0\t0/1\n",
        );
        // bgzip + tabix index
        let status = Command::new("bgzip").arg(&vcf).status().expect("bgzip");
        assert!(status.success(), "bgzip failed");
        let vcf_gz = vcf.with_extension("vcf.gz");
        let status = Command::new("tabix")
            .arg("-p")
            .arg("vcf")
            .arg(&vcf_gz)
            .status()
            .expect("tabix");
        assert!(status.success(), "tabix failed");

        let store = VcfStore::load(&vcf_gz).unwrap();

        // Several regions, 1-based inclusive. bcftools view -r default overlap=1.
        for (start1, end1) in [(1i64, 15), (6, 20), (10, 30), (1, 100), (26, 29)] {
            let mine = store.query("chr1", start1 - 1, end1 - 1, 1).len();

            let out = Command::new("bcftools")
                .args(["view", "-r", &format!("chr1:{}-{}", start1, end1)])
                .arg(&vcf_gz)
                .output()
                .expect("bcftools");
            assert!(
                out.status.success(),
                "bcftools view failed: {}",
                String::from_utf8_lossy(&out.stderr)
            );
            let theirs = String::from_utf8_lossy(&out.stdout)
                .lines()
                .filter(|l| !l.is_empty() && !l.starts_with('#'))
                .count();
            assert_eq!(
                mine, theirs,
                "region chr1:{}-{}: mine={} bcftools={}",
                start1, end1, mine, theirs
            );
        }
    }
}
