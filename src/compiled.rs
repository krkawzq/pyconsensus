//! compiled — preclassified VCF record metadata for fastpath-first execution.
//!
//! This is the first migration step from interpretive `VcfRecord` execution to
//! a compiled VCF representation. It intentionally stays compatible with the
//! existing `VcfStore`: every loaded record receives a compact classification
//! and per-allele operation table, while the legacy fields remain available for
//! fallback correctness.

use smallvec::SmallVec;

#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RecordKind {
    RefOnly = 0,
    Snp1 = 1,
    SameLen = 2,
    NormInsertion = 3,
    NormDeletion = 4,
    SimpleIndel = 5,
    SymbolicDel = 6,
    GvcfBlock = 7,
    Complex = 8,
}

impl RecordKind {
    pub const COUNT: usize = 9;

    #[inline]
    pub fn as_usize(self) -> usize {
        self as usize
    }

    pub fn name(self) -> &'static str {
        match self {
            RecordKind::RefOnly => "RefOnly",
            RecordKind::Snp1 => "Snp1",
            RecordKind::SameLen => "SameLen",
            RecordKind::NormInsertion => "NormInsertion",
            RecordKind::NormDeletion => "NormDeletion",
            RecordKind::SimpleIndel => "SimpleIndel",
            RecordKind::SymbolicDel => "SymbolicDel",
            RecordKind::GvcfBlock => "GvcfBlock",
            RecordKind::Complex => "Complex",
        }
    }
}

#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AlleleOpKind {
    Ref = 0,
    SameLen = 1,
    Insert = 2,
    Delete = 3,
    Replace = 4,
    SymbolicDel = 5,
    GvcfRefBlock = 6,
    Missing = 7,
    Unsupported = 8,
}

impl AlleleOpKind {
    pub const COUNT: usize = 9;

    #[inline]
    pub fn as_usize(self) -> usize {
        self as usize
    }

    pub fn name(self) -> &'static str {
        match self {
            AlleleOpKind::Ref => "Ref",
            AlleleOpKind::SameLen => "SameLen",
            AlleleOpKind::Insert => "Insert",
            AlleleOpKind::Delete => "Delete",
            AlleleOpKind::Replace => "Replace",
            AlleleOpKind::SymbolicDel => "SymbolicDel",
            AlleleOpKind::GvcfRefBlock => "GvcfRefBlock",
            AlleleOpKind::Missing => "Missing",
            AlleleOpKind::Unsupported => "Unsupported",
        }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RecordFlags {
    bits: u16,
}

impl RecordFlags {
    pub const BIALLELIC: u16 = 1 << 0;
    pub const MULTI_ALLELIC: u16 = 1 << 1;
    pub const HAS_SYMBOLIC: u16 = 1 << 2;
    pub const HAS_STAR: u16 = 1 << 3;
    pub const HAS_LEN_CHANGE: u16 = 1 << 4;
    pub const ALL_ALT_SAME_LEN: u16 = 1 << 5;
    pub const ALL_ALT_FASTPATH_ELIGIBLE: u16 = 1 << 6;

    #[inline]
    pub fn from_bits(bits: u16) -> Self {
        RecordFlags { bits }
    }

    #[inline]
    pub fn set(&mut self, bit: u16) {
        self.bits |= bit;
    }

    #[inline]
    pub fn contains(self, bit: u16) -> bool {
        self.bits & bit != 0
    }

    #[inline]
    pub fn bits(self) -> u16 {
        self.bits
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AlleleOp {
    pub kind: AlleleOpKind,
    pub ref_len: u32,
    pub alt_len: u32,
    pub trim_beg: u16,
    pub len_diff: i32,
    pub case_flags: u8,
}

pub const ALLELE_HAS_ASCII_LOWER: u8 = 1 << 0;
pub const ALLELE_HAS_ASCII_UPPER: u8 = 1 << 1;

impl AlleleOp {
    pub fn ref_op(ref_len: u32, ref_allele: &[u8]) -> Self {
        AlleleOp {
            kind: AlleleOpKind::Ref,
            ref_len,
            alt_len: ref_len,
            trim_beg: 0,
            len_diff: 0,
            case_flags: allele_case_flags(ref_allele),
        }
    }

    #[inline]
    pub fn is_same_len_fastpath(&self) -> bool {
        matches!(self.kind, AlleleOpKind::Ref | AlleleOpKind::SameLen)
            && self.len_diff == 0
            && self.ref_len == self.alt_len
    }

    #[inline]
    pub fn is_edit_script_fastpath(&self) -> bool {
        matches!(
            self.kind,
            AlleleOpKind::SameLen | AlleleOpKind::Insert | AlleleOpKind::Delete
        )
    }

    #[inline]
    pub fn has_ascii_lowercase(&self) -> bool {
        self.case_flags & ALLELE_HAS_ASCII_LOWER != 0
    }

    #[inline]
    pub fn has_ascii_uppercase(&self) -> bool {
        self.case_flags & ALLELE_HAS_ASCII_UPPER != 0
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CompiledRecord {
    pub kind: RecordKind,
    pub flags: RecordFlags,
    pub ops: SmallVec<[AlleleOp; 2]>,
}

impl CompiledRecord {
    pub fn from_alleles(rlen: i32, alleles: &[SmallVec<[u8; 16]>]) -> Self {
        let ref_len = rlen.max(0) as u32;
        let ref_allele = alleles.first().map(|a| &a[..]).unwrap_or(&[]);
        let mut flags = RecordFlags::default();
        if alleles.len() == 2 {
            flags.set(RecordFlags::BIALLELIC);
        } else if alleles.len() > 2 {
            flags.set(RecordFlags::MULTI_ALLELIC);
        }

        let mut ops: SmallVec<[AlleleOp; 2]> = SmallVec::with_capacity(alleles.len().max(1));
        ops.push(AlleleOp::ref_op(ref_len, ref_allele));
        for alt in alleles.iter().skip(1) {
            let op = compile_alt_op(ref_allele, ref_len, alt);
            if matches!(
                op.kind,
                AlleleOpKind::SymbolicDel | AlleleOpKind::GvcfRefBlock
            ) {
                flags.set(RecordFlags::HAS_SYMBOLIC);
            }
            if alt.as_slice() == b"*" {
                flags.set(RecordFlags::HAS_STAR);
            }
            if op.len_diff != 0 {
                flags.set(RecordFlags::HAS_LEN_CHANGE);
            }
            ops.push(op);
        }

        let n_alt = ops.len().saturating_sub(1);
        if n_alt > 0 && ops.iter().skip(1).all(AlleleOp::is_same_len_fastpath) {
            flags.set(RecordFlags::ALL_ALT_SAME_LEN);
            flags.set(RecordFlags::ALL_ALT_FASTPATH_ELIGIBLE);
        } else if n_alt > 0 && ops.iter().skip(1).all(AlleleOp::is_edit_script_fastpath) {
            flags.set(RecordFlags::ALL_ALT_FASTPATH_ELIGIBLE);
        }

        let kind = classify_record_kind(&ops, n_alt);
        CompiledRecord { kind, flags, ops }
    }

    #[inline]
    pub fn allele_op(&self, ialt: usize) -> Option<&AlleleOp> {
        self.ops.get(ialt)
    }

    #[inline]
    pub fn same_len_allele(&self, ialt: usize) -> bool {
        self.allele_op(ialt)
            .map(AlleleOp::is_same_len_fastpath)
            .unwrap_or(false)
    }
}

impl Default for CompiledRecord {
    fn default() -> Self {
        CompiledRecord {
            kind: RecordKind::Complex,
            flags: RecordFlags::default(),
            ops: SmallVec::new(),
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct VcfCompileStats {
    pub records_total: u64,
    pub kind_counts: [u64; RecordKind::COUNT],
    pub allele_op_counts: [u64; AlleleOpKind::COUNT],
    pub biallelic_records: u64,
    pub multi_allelic_records: u64,
    pub gt_records: u64,
    pub biallelic_phased_diploid_records: u64,
    pub compact_gt_records: u64,
    pub biallelic_gt_bitset_records: u64,
    pub missing_gt_records: u64,
}

impl VcfCompileStats {
    pub fn observe_record(&mut self, record: &CompiledRecord) {
        self.records_total += 1;
        self.kind_counts[record.kind.as_usize()] += 1;
        for op in &record.ops {
            self.allele_op_counts[op.kind.as_usize()] += 1;
        }
        if record.flags.contains(RecordFlags::BIALLELIC) {
            self.biallelic_records += 1;
        }
        if record.flags.contains(RecordFlags::MULTI_ALLELIC) {
            self.multi_allelic_records += 1;
        }
    }

    pub fn observe_gt(
        &mut self,
        has_gt: bool,
        is_biallelic_phased_diploid: bool,
        has_compact_gt: bool,
        has_biallelic_gt_bitset: bool,
        has_missing: bool,
    ) {
        if has_gt {
            self.gt_records += 1;
        }
        if is_biallelic_phased_diploid {
            self.biallelic_phased_diploid_records += 1;
        }
        if has_compact_gt {
            self.compact_gt_records += 1;
        }
        if has_biallelic_gt_bitset {
            self.biallelic_gt_bitset_records += 1;
        }
        if has_missing {
            self.missing_gt_records += 1;
        }
    }

    pub fn kind_count(&self, kind: RecordKind) -> u64 {
        self.kind_counts[kind.as_usize()]
    }

    pub fn allele_op_count(&self, kind: AlleleOpKind) -> u64 {
        self.allele_op_counts[kind.as_usize()]
    }

    pub fn summary_lines(&self) -> Vec<String> {
        let mut lines = vec![format!("records_total={}", self.records_total)];
        for kind in [
            RecordKind::RefOnly,
            RecordKind::Snp1,
            RecordKind::SameLen,
            RecordKind::NormInsertion,
            RecordKind::NormDeletion,
            RecordKind::SimpleIndel,
            RecordKind::SymbolicDel,
            RecordKind::GvcfBlock,
            RecordKind::Complex,
        ] {
            lines.push(format!(
                "record_kind.{}={}",
                kind.name(),
                self.kind_count(kind)
            ));
        }
        lines.push(format!("biallelic_records={}", self.biallelic_records));
        lines.push(format!(
            "multi_allelic_records={}",
            self.multi_allelic_records
        ));
        lines.push(format!("gt_records={}", self.gt_records));
        lines.push(format!(
            "biallelic_phased_diploid_records={}",
            self.biallelic_phased_diploid_records
        ));
        lines.push(format!("compact_gt_records={}", self.compact_gt_records));
        lines.push(format!(
            "biallelic_gt_bitset_records={}",
            self.biallelic_gt_bitset_records
        ));
        lines.push(format!("missing_gt_records={}", self.missing_gt_records));
        lines
    }
}

fn compile_alt_op(ref_allele: &[u8], ref_len: u32, alt: &[u8]) -> AlleleOp {
    let alt_len = alt.len() as u32;
    let case_flags = allele_case_flags(alt);
    if alt.eq_ignore_ascii_case(b"<DEL>") {
        return AlleleOp {
            kind: AlleleOpKind::SymbolicDel,
            ref_len,
            alt_len: 1,
            trim_beg: 1,
            len_diff: 1 - ref_len as i32,
            case_flags,
        };
    }
    if alt.eq_ignore_ascii_case(b"<*>") || alt.eq_ignore_ascii_case(b"<NON_REF>") {
        return AlleleOp {
            kind: AlleleOpKind::GvcfRefBlock,
            ref_len,
            alt_len,
            trim_beg: 0,
            len_diff: 0,
            case_flags,
        };
    }
    if alt.starts_with(b"<") || alt == b"*" || alt.is_empty() {
        return AlleleOp {
            kind: AlleleOpKind::Unsupported,
            ref_len,
            alt_len,
            trim_beg: 0,
            len_diff: alt_len as i32 - ref_len as i32,
            case_flags,
        };
    }

    let trim_beg = if !ref_allele.is_empty() && ref_allele.first() == alt.first() {
        1
    } else {
        0
    };
    let len_diff = alt_len as i32 - ref_len as i32;
    let kind = if alt_len == ref_len {
        AlleleOpKind::SameLen
    } else if trim_beg == 1 && alt_len > ref_len {
        AlleleOpKind::Insert
    } else if trim_beg == 1 && alt_len < ref_len {
        AlleleOpKind::Delete
    } else if alt_len != ref_len {
        AlleleOpKind::Replace
    } else {
        AlleleOpKind::Unsupported
    };

    AlleleOp {
        kind,
        ref_len,
        alt_len,
        trim_beg,
        len_diff,
        case_flags,
    }
}

pub fn allele_case_flags(bytes: &[u8]) -> u8 {
    let mut flags = 0u8;
    for &b in bytes {
        if b.is_ascii_lowercase() {
            flags |= ALLELE_HAS_ASCII_LOWER;
        } else if b.is_ascii_uppercase() {
            flags |= ALLELE_HAS_ASCII_UPPER;
        }
    }
    flags
}

fn classify_record_kind(ops: &[AlleleOp], n_alt: usize) -> RecordKind {
    if n_alt == 0 {
        return RecordKind::RefOnly;
    }
    let alts = &ops[1..];
    if alts
        .iter()
        .all(|op| op.kind == AlleleOpKind::SameLen && op.ref_len == 1)
    {
        return RecordKind::Snp1;
    }
    if alts.iter().all(|op| op.kind == AlleleOpKind::SameLen) {
        return RecordKind::SameLen;
    }
    if alts.iter().all(|op| op.kind == AlleleOpKind::Insert) {
        return RecordKind::NormInsertion;
    }
    if alts.iter().all(|op| op.kind == AlleleOpKind::Delete) {
        return RecordKind::NormDeletion;
    }
    if alts.iter().all(|op| op.kind == AlleleOpKind::SymbolicDel) {
        return RecordKind::SymbolicDel;
    }
    if alts.iter().all(|op| op.kind == AlleleOpKind::GvcfRefBlock) {
        return RecordKind::GvcfBlock;
    }
    if alts.iter().all(|op| {
        matches!(
            op.kind,
            AlleleOpKind::Insert | AlleleOpKind::Delete | AlleleOpKind::Replace
        )
    }) {
        return RecordKind::SimpleIndel;
    }
    RecordKind::Complex
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alleles(xs: &[&[u8]]) -> Vec<SmallVec<[u8; 16]>> {
        xs.iter().map(|x| SmallVec::from_slice(x)).collect()
    }

    #[test]
    fn classifies_common_record_shapes() {
        let snp = CompiledRecord::from_alleles(1, &alleles(&[b"A", b"G"]));
        assert_eq!(snp.kind, RecordKind::Snp1);
        assert!(snp.same_len_allele(1));

        let mnp = CompiledRecord::from_alleles(2, &alleles(&[b"AC", b"GT"]));
        assert_eq!(mnp.kind, RecordKind::SameLen);
        assert!(mnp.flags.contains(RecordFlags::ALL_ALT_SAME_LEN));

        let ins = CompiledRecord::from_alleles(1, &alleles(&[b"A", b"AT"]));
        assert_eq!(ins.kind, RecordKind::NormInsertion);
        assert_eq!(ins.ops[1].len_diff, 1);

        let del = CompiledRecord::from_alleles(2, &alleles(&[b"AC", b"A"]));
        assert_eq!(del.kind, RecordKind::NormDeletion);
        assert_eq!(del.ops[1].len_diff, -1);

        let sym = CompiledRecord::from_alleles(2, &alleles(&[b"AC", b"<DEL>"]));
        assert_eq!(sym.kind, RecordKind::SymbolicDel);

        let gvcf = CompiledRecord::from_alleles(1, &alleles(&[b"A", b"<NON_REF>"]));
        assert_eq!(gvcf.kind, RecordKind::GvcfBlock);
    }
}
