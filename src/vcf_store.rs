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

use crate::compiled::{CompiledRecord, VcfCompileStats};
use crate::htslib_ffi as ffi;
use crate::planner::{plan_region, PlanOptions, RegionPlan};
use smallvec::SmallVec;
use std::collections::HashMap;
use std::ffi::CString;
use std::fs::{self, File};
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::os::raw::{c_int, c_void};
use std::path::{Path, PathBuf};
use std::time::UNIX_EPOCH;

const CVCF_MAGIC: &[u8; 8] = b"CVCF0001";
const CVCF_VERSION: u32 = 2;

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
    /// Per-sample GT; `len == n_sample` when the VCF has GT, else empty.
    /// Ploidy is variable per sample.
    pub gt: Vec<SmallVec<[GtAllele; 2]>>,
    /// Bitset GT fastpath for biallelic phased diploid records.
    pub gt_bits: Option<BiallelicPhasedGtBits>,
    /// `bcf_get_variant_types` bitmask (VCF_SNP|MNP|INDEL|...), precomputed.
    pub var_type: i32,
    /// Preclassified record/allele metadata used by fastpath dispatch.
    pub compiled: CompiledRecord,
}

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
        let rid = match self.seq_names.get(chr) {
            Some(r) => *r,
            None => return Vec::new(),
        };
        let idx = match self.by_rid.get(&rid) {
            Some(v) => v,
            None => return Vec::new(),
        };
        if idx.is_empty() {
            return Vec::new();
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
        let base_cap = hi.saturating_sub(lo_pos);
        let span_cap = lo_pos.saturating_sub(first_spanning);
        let mut out: Vec<&VcfRecord> = Vec::with_capacity(base_cap + span_cap);
        match overlap {
            0 => {
                for k in lo_pos..hi {
                    out.push(&self.records[idx[k] as usize]);
                }
            }
            1 | 2 => {
                // Records with pos < start that still span into [start, end].
                if lo_pos > 0 {
                    for k in first_spanning..lo_pos {
                        let rec = &self.records[idx[k] as usize];
                        if rec.pos <= end && rec.ref_end() >= start {
                            out.push(rec);
                        }
                    }
                }
                // Records with pos in [start, end]: ref_end >= pos >= start, always overlap.
                for k in lo_pos..hi {
                    out.push(&self.records[idx[k] as usize]);
                }
            }
            _ => return Vec::new(),
        }
        out
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
        for record_idx in 0..n_records {
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

            let compiled = CompiledRecord::from_alleles(rlen, &alleles);
            store.compile_stats.observe_record(&compiled);
            let (has_gt, is_biallelic_phased_diploid, has_missing_gt) =
                gt_compile_stats(n_allele, &gt);
            let gt_bits = BiallelicPhasedGtBits::from_gt(n_allele, store.n_sample as usize, &gt);
            store.compile_stats.observe_gt(
                has_gt,
                is_biallelic_phased_diploid,
                gt_bits.is_some(),
                has_missing_gt,
            );
            store.has_gt |= has_gt;

            store.by_rid.entry(rid).or_default().push(record_idx as u32);
            store.records.push(VcfRecord {
                pos,
                rlen,
                rid,
                alleles,
                gt,
                gt_bits,
                var_type,
                compiled,
            });
        }
        store.rebuild_pmax_end();
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
        }
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
                let mut gt: Vec<SmallVec<[GtAllele; 2]>> = Vec::new();
                let mut has_missing_gt = false;
                let mut is_biallelic_phased_diploid = n_allele == 2;
                if ngt > 0 && self.n_sample > 0 && !gt_buf.is_null() {
                    self.has_gt = true;
                    gt.reserve_exact(self.n_sample as usize);
                    // bcf_get_format_values returns nsmpl*max_ploidy (the total
                    // number of int32 values), so per-sample ploidy is:
                    let ploidy = (ngt as usize) / (self.n_sample as usize);
                    let base = gt_buf as *const i32;
                    for s in 0..self.n_sample as usize {
                        let mut alleles_g: SmallVec<[GtAllele; 2]> = SmallVec::new();
                        for j in 0..ploidy {
                            let raw = *base.add(s * ploidy + j);
                            if raw == ffi::BCF_INT32_VECTOR_END {
                                break;
                            }
                            // Note: BCF_INT32_MISSING (-2147483648) should not
                            // appear in GT-encoded arrays (missing there = 0),
                            // but guard anyway.
                            let allele = if ffi::gt_is_missing(raw) {
                                has_missing_gt = true;
                                None
                            } else {
                                Some(ffi::gt_allele(raw))
                            };
                            alleles_g.push(GtAllele {
                                allele,
                                phased: ffi::gt_is_phased(raw),
                                raw,
                            });
                        }
                        if alleles_g.len() != 2
                            || alleles_g.iter().any(|a| a.allele.is_none())
                            || !alleles_g.iter().all(|a| a.phased)
                        {
                            is_biallelic_phased_diploid = false;
                        }
                        gt.push(alleles_g);
                    }
                } else {
                    is_biallelic_phased_diploid = false;
                }

                let var_type = ffi::bcf_get_variant_types(rec);
                let compiled = CompiledRecord::from_alleles(rlen, &alleles);
                let gt_bits =
                    BiallelicPhasedGtBits::from_gt(n_allele as usize, self.n_sample as usize, &gt);
                self.compile_stats.observe_record(&compiled);
                self.compile_stats.observe_gt(
                    !gt.is_empty(),
                    is_biallelic_phased_diploid,
                    gt_bits.is_some(),
                    has_missing_gt,
                );

                let rid_bucket = self.by_rid.entry(rid).or_default();
                rid_bucket.push(self.records.len() as u32);
                self.records.push(VcfRecord {
                    pos,
                    rlen,
                    rid,
                    alleles,
                    gt,
                    gt_bits,
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

fn read_u8<R: Read>(r: &mut R) -> io::Result<u8> {
    let mut b = [0u8; 1];
    r.read_exact(&mut b)?;
    Ok(b[0])
}

fn write_u8<W: Write>(w: &mut W, v: u8) -> io::Result<()> {
    w.write_all(&[v])
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

        let r0 = &store.records()[0];
        assert_eq!(r0.pos, 9); // 1-based 10 -> 0-based 9
        assert_eq!(r0.rlen, 1);
        assert_eq!(r0.alleles.len(), 2);
        assert_eq!(&r0.alleles[0][..], b"G");
        assert_eq!(&r0.alleles[1][..], b"A");
        assert_eq!(r0.gt.len(), 2);
        // S1 = 0|1 (phased): allele0=REF(0) phased, allele1=ALT(1) phased
        assert_eq!(r0.gt[0][0].allele, Some(0));
        assert!(r0.gt[0][0].phased);
        assert_eq!(r0.gt[0][1].allele, Some(1));
        // S2 = 1/1 (unphased): both ALT, unphased
        assert_eq!(r0.gt[1][0].allele, Some(1));
        assert!(!r0.gt[1][0].phased);

        // missing GT ./.
        let r1 = &store.records()[1];
        assert!(r1.gt[1][0].allele.is_none());
    }

    #[test]
    fn writes_and_reads_owned_cvcf_cache() {
        let vcf = write_vcf(
            "cache_roundtrip",
            "chr1\t10\t.\tG\tA\t.\t.\t.\tGT\t0|1\t1/1\n\
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
        let q = cached.query("chr1", 0, 25, 1);
        assert_eq!(q.len(), 2);
        assert_eq!(&q[0].alleles[1][..], b"A");
        assert!(q[1].gt[1][0].allele.is_none());
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
