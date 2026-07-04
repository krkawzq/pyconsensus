//! planner — region-level fastpath classification.
//!
//! A `RegionPlan` is deliberately conservative: it never proves correctness by
//! itself, it only identifies the fastest lane that is safe to try before
//! falling back to the legacy state machine. Allele selection is still runtime
//! dependent, so this planner classifies by record/allele capabilities rather
//! than by a concrete sample's final alleles.

use crate::compiled::{RecordFlags, RecordKind};
use crate::stats::{FallbackReason, FastPathLane};
use crate::vcf_store::{RecordSet, VcfRecord};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PlanOptions {
    pub chain: bool,
    pub mark_del: bool,
    pub mark_ins: bool,
    pub mark_snv: bool,
    pub mask: bool,
    pub mask_skips_variants: bool,
    pub mask_overlaps_variant: bool,
    pub absent: bool,
}

impl PlanOptions {
    #[inline]
    pub fn has_mark(self) -> bool {
        self.mark_del || self.mark_ins || self.mark_snv
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RegionPlan {
    pub lane: FastPathLane,
    pub records_total: usize,
    pub same_len_records: usize,
    pub edit_script_records: usize,
    pub fallback_records: usize,
    pub fallback_reasons: Vec<FallbackReason>,
    pub mask_overlap_known: bool,
    pub mask_overlaps_variant: bool,
}

impl RegionPlan {
    pub fn empty() -> Self {
        RegionPlan {
            lane: FastPathLane::EmptyRegion,
            records_total: 0,
            same_len_records: 0,
            edit_script_records: 0,
            fallback_records: 0,
            fallback_reasons: Vec::new(),
            mask_overlap_known: false,
            mask_overlaps_variant: false,
        }
    }

    pub fn needs_fallback(reason: FallbackReason, records_total: usize) -> Self {
        let mask_overlaps_variant = reason == FallbackReason::MaskOverlap;
        RegionPlan {
            lane: FastPathLane::FallbackStateMachine,
            records_total,
            same_len_records: 0,
            edit_script_records: 0,
            fallback_records: records_total,
            fallback_reasons: vec![reason],
            mask_overlap_known: mask_overlaps_variant,
            mask_overlaps_variant,
        }
    }
}

pub fn plan_region(records: &[&VcfRecord], opts: PlanOptions) -> RegionPlan {
    plan_region_set(&RecordSet::from_ref_slice(records), opts)
}

pub fn plan_region_set(records: &RecordSet<'_>, opts: PlanOptions) -> RegionPlan {
    if records.is_empty() {
        return RegionPlan::empty();
    }

    if opts.mask && opts.mask_skips_variants && opts.mask_overlaps_variant {
        return RegionPlan::needs_fallback(FallbackReason::MaskOverlap, records.len());
    }
    let mask_overlap_known = opts.mask && opts.mask_skips_variants;
    let mut same_len_records = 0usize;
    let mut edit_script_records = 0usize;
    let mut fallback_records = 0usize;
    let mut fallback_reasons = Vec::new();

    for meta in records.iter_meta() {
        if meta.flags.contains(RecordFlags::HAS_SYMBOLIC) && meta.kind != RecordKind::SymbolicDel {
            fallback_records += 1;
            push_unique(&mut fallback_reasons, FallbackReason::SymbolicAllele);
            continue;
        }
        if meta.flags.contains(RecordFlags::HAS_STAR) {
            fallback_records += 1;
            push_unique(&mut fallback_reasons, FallbackReason::ComplexAllele);
            continue;
        }

        if meta.kind == RecordKind::RefOnly || meta.flags.contains(RecordFlags::ALL_ALT_SAME_LEN) {
            same_len_records += 1;
        } else if meta.flags.contains(RecordFlags::ALL_ALT_FASTPATH_ELIGIBLE) {
            edit_script_records += 1;
        } else {
            fallback_records += 1;
            push_unique(&mut fallback_reasons, FallbackReason::ComplexAllele);
        }
    }

    let lane = if fallback_records > 0 {
        FastPathLane::FallbackStateMachine
    } else if edit_script_records == 0 {
        FastPathLane::SameLenOnly
    } else if same_len_records == 0 {
        FastPathLane::NormalizedEditScript
    } else {
        FastPathLane::MixedSimpleEdits
    };

    RegionPlan {
        lane,
        records_total: records.len(),
        same_len_records,
        edit_script_records,
        fallback_records,
        fallback_reasons,
        mask_overlap_known,
        mask_overlaps_variant: opts.mask_overlaps_variant,
    }
}

fn push_unique(xs: &mut Vec<FallbackReason>, x: FallbackReason) {
    if !xs.contains(&x) {
        xs.push(x);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compiled::CompiledRecord;
    use crate::vcf_store::VcfRecord;
    use smallvec::SmallVec;

    fn rec(pos: i64, rlen: i32, alleles: &[&[u8]]) -> VcfRecord {
        let alleles: Vec<SmallVec<[u8; 16]>> =
            alleles.iter().map(|a| SmallVec::from_slice(a)).collect();
        let compiled = CompiledRecord::from_alleles(rlen, &alleles);
        VcfRecord {
            pos,
            rlen,
            rid: 0,
            alleles,
            gt: Vec::new(),
            gt_compact: None,
            gt_bits: None,
            var_type: 0,
            compiled,
        }
    }

    #[test]
    fn plans_same_len_and_edit_script_regions() {
        let snp = rec(1, 1, &[b"A", b"G"]);
        let mnp = rec(3, 2, &[b"AC", b"GT"]);
        let plan = plan_region(&[&snp, &mnp], PlanOptions::default());
        assert_eq!(plan.lane, FastPathLane::SameLenOnly);
        assert_eq!(plan.same_len_records, 2);

        let ins = rec(5, 1, &[b"A", b"AT"]);
        let del = rec(8, 2, &[b"AC", b"A"]);
        let plan = plan_region(&[&ins, &del], PlanOptions::default());
        assert_eq!(plan.lane, FastPathLane::NormalizedEditScript);
        assert_eq!(plan.edit_script_records, 2);

        let marked_plan = plan_region(
            &[&ins, &del],
            PlanOptions {
                mark_del: true,
                mark_ins: true,
                mark_snv: true,
                ..Default::default()
            },
        );
        assert_eq!(marked_plan.lane, FastPathLane::NormalizedEditScript);
        assert_eq!(marked_plan.edit_script_records, 2);

        let plan = plan_region(&[&snp, &ins], PlanOptions::default());
        assert_eq!(plan.lane, FastPathLane::MixedSimpleEdits);
    }

    #[test]
    fn plans_symbolic_del_as_edit_script_and_chain_keeps_simple_lanes() {
        let sym = rec(1, 2, &[b"AC", b"<DEL>"]);
        let plan = plan_region(&[&sym], PlanOptions::default());
        assert_eq!(plan.lane, FastPathLane::NormalizedEditScript);
        assert_eq!(plan.edit_script_records, 1);

        let gvcf = rec(1, 2, &[b"AC", b"<NON_REF>"]);
        let plan = plan_region(&[&gvcf], PlanOptions::default());
        assert_eq!(plan.lane, FastPathLane::FallbackStateMachine);
        assert!(plan
            .fallback_reasons
            .contains(&FallbackReason::SymbolicAllele));

        let snp = rec(1, 1, &[b"A", b"G"]);
        let plan = plan_region(
            &[&snp],
            PlanOptions {
                chain: true,
                ..Default::default()
            },
        );
        assert_eq!(plan.lane, FastPathLane::SameLenOnly);

        let ins = rec(5, 1, &[b"A", b"AT"]);
        let del = rec(8, 2, &[b"AC", b"A"]);
        let plan = plan_region(
            &[&ins, &del],
            PlanOptions {
                chain: true,
                ..Default::default()
            },
        );
        assert_eq!(plan.lane, FastPathLane::NormalizedEditScript);
    }

    #[test]
    fn records_mask_overlap_metadata_in_plan() {
        let snp = rec(1, 1, &[b"A", b"G"]);
        let no_overlap = plan_region(
            &[&snp],
            PlanOptions {
                mask: true,
                mask_skips_variants: true,
                mask_overlaps_variant: false,
                ..Default::default()
            },
        );
        assert_eq!(no_overlap.lane, FastPathLane::SameLenOnly);
        assert!(no_overlap.mask_overlap_known);
        assert!(!no_overlap.mask_overlaps_variant);

        let overlap = plan_region(
            &[&snp],
            PlanOptions {
                mask: true,
                mask_skips_variants: true,
                mask_overlaps_variant: true,
                ..Default::default()
            },
        );
        assert_eq!(overlap.lane, FastPathLane::FallbackStateMachine);
        assert!(overlap.mask_overlap_known);
        assert!(overlap.mask_overlaps_variant);
        assert!(overlap
            .fallback_reasons
            .contains(&FallbackReason::MaskOverlap));
    }
}
