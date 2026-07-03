//! planner — region-level fastpath classification.
//!
//! A `RegionPlan` is deliberately conservative: it never proves correctness by
//! itself, it only identifies the fastest lane that is safe to try before
//! falling back to the legacy state machine. Allele selection is still runtime
//! dependent, so this planner classifies by record/allele capabilities rather
//! than by a concrete sample's final alleles.

use crate::compiled::{AlleleOp, AlleleOpKind, RecordFlags};
use crate::stats::{FallbackReason, FastPathLane};
use crate::vcf_store::VcfRecord;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PlanOptions {
    pub chain: bool,
    pub mark_del: bool,
    pub mark_ins: bool,
    pub mark_snv: bool,
    pub mask: bool,
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
        }
    }

    pub fn needs_fallback(reason: FallbackReason, records_total: usize) -> Self {
        RegionPlan {
            lane: FastPathLane::FallbackStateMachine,
            records_total,
            same_len_records: 0,
            edit_script_records: 0,
            fallback_records: records_total,
            fallback_reasons: vec![reason],
        }
    }
}

pub fn plan_region(records: &[&VcfRecord], opts: PlanOptions) -> RegionPlan {
    if records.is_empty() {
        return RegionPlan::empty();
    }

    if opts.chain {
        return RegionPlan::needs_fallback(FallbackReason::ChainEnabled, records.len());
    }
    if opts.mask {
        return RegionPlan::needs_fallback(FallbackReason::MaskOverlap, records.len());
    }
    if opts.has_mark()
        && records
            .iter()
            .any(|r| r.compiled.flags.contains(RecordFlags::HAS_LEN_CHANGE))
    {
        return RegionPlan::needs_fallback(FallbackReason::MarkEnabled, records.len());
    }

    let mut same_len_records = 0usize;
    let mut edit_script_records = 0usize;
    let mut fallback_records = 0usize;
    let mut fallback_reasons = Vec::new();

    for rec in records {
        if rec.compiled.flags.contains(RecordFlags::HAS_SYMBOLIC) {
            fallback_records += 1;
            push_unique(&mut fallback_reasons, FallbackReason::SymbolicAllele);
            continue;
        }
        if rec.compiled.flags.contains(RecordFlags::HAS_STAR) {
            fallback_records += 1;
            push_unique(&mut fallback_reasons, FallbackReason::ComplexAllele);
            continue;
        }

        let alts = rec.compiled.ops.iter().skip(1);
        let mut n_alt = 0usize;
        let mut all_same_len = true;
        let mut all_edit_script = true;
        for op in alts {
            n_alt += 1;
            if !op.is_same_len_fastpath() {
                all_same_len = false;
            }
            if !is_normalized_edit_script_op(op) {
                all_edit_script = false;
            }
        }

        if n_alt == 0 || all_same_len {
            same_len_records += 1;
        } else if all_edit_script {
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
    }
}

fn is_normalized_edit_script_op(op: &AlleleOp) -> bool {
    matches!(
        op.kind,
        AlleleOpKind::SameLen | AlleleOpKind::Insert | AlleleOpKind::Delete
    )
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

        let plan = plan_region(&[&snp, &ins], PlanOptions::default());
        assert_eq!(plan.lane, FastPathLane::MixedSimpleEdits);
    }

    #[test]
    fn plans_symbolic_and_chain_as_fallback() {
        let sym = rec(1, 2, &[b"AC", b"<DEL>"]);
        let plan = plan_region(&[&sym], PlanOptions::default());
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
        assert_eq!(plan.lane, FastPathLane::FallbackStateMachine);
        assert!(plan
            .fallback_reasons
            .contains(&FallbackReason::ChainEnabled));
    }
}
