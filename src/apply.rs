//! apply — consensus apply state machine (Rust rewrite of bcftools consensus).
//!
//! Faithful port of `apply_variant()` (bcftools/consensus.c:583) and
//! helpers `apply_absent` / `freeze_ref`. Line numbers in comments refer to
//! that file.
//!
//! Key difference from bcftools (docs/implementation_plan.md §0.3): bcftools
//! streams a `fa_buf` across a whole chromosome (flushing as it goes); we load
//! the entire region into `state.buf` once and return it whole. The
//! `ApplyState` machine (`mod_off`/`frz_pos`/`prev_base`/`frz_mod`/...) is
//! replicated exactly — only the streaming `flush_fa_buffer` is absent.
//!
//! M2 scope: no sample, no `-I` → `ialt = 1` (apply the first ALT; for
//! biallelic data this is the only ALT). Mask / chain / mark / absent /
//! missing hooks are present but inert until M4 wires real options.

use crate::chain::Chain;
use crate::compiled::AlleleOpKind;
use crate::haplotype::{select_allele, AlleleSelection, SampleMode};
use crate::mask::Mask;
use crate::planner::{plan_region, PlanOptions, RegionPlan};
use crate::stats::{FallbackReason, FastPathLane, RuntimeStats};
use crate::vcf_store::VcfRecord;
use smallvec::SmallVec;
use std::rc::Rc;

type AlleleBuf = SmallVec<[u8; 64]>;

pub const TO_UPPER: i8 = 1;
pub const TO_LOWER: i8 = 2;

/// Options for a consensus apply pass. Maps 1:1 to bcftools consensus flags.
#[derive(Default)]
pub struct ApplyOptions {
    /// `-a` absent char: fill positions in the region not covered by any VCF record.
    pub absent_allele: Option<u8>,
    /// `-M` missing char: emit this for missing GT (`.`) instead of skipping.
    pub missing_allele: Option<u8>,
    /// `--mark-del` char. None = delete bases normally.
    pub mark_del: Option<u8>,
    /// `--mark-ins`: Some(TO_UPPER|TO_LOWER) or Some(char).
    pub mark_ins: Option<u8>,
    /// `--mark-snv`: Some(TO_UPPER|TO_LOWER) or Some(char).
    pub mark_snv: Option<u8>,
    /// How to pick the allele per record. Defaults to ApplyAllAlt (no `-s`, no `-I`).
    pub sample_mode: SampleMode,
    /// `-m` mask (char-mode skips overlapping variants). Rc lets a region group
    /// reuse one mask instance sequentially without sharing htslib iterator
    /// state across threads.
    pub mask: Option<Rc<Mask>>,
}

impl ApplyOptions {
    /// The M2 default: no sample, no `-I` → apply the first ALT.
    pub fn apply_all_alt() -> Self {
        ApplyOptions {
            sample_mode: SampleMode::ApplyAllAlt,
            ..Default::default()
        }
    }
}

/// The apply state machine. `buf` is the in-progress consensus (starts as the
/// fetched ref subsequence, plus-strand).
pub struct ApplyState {
    /// region start (0-based), fixed — `fa_ori_pos`.
    pub ori_pos: i64,
    /// modified-vs-ori offset — `fa_mod_off`.
    pub mod_off: i64,
    /// frozen original coord (last applied variant end) — `fa_frz_pos`. -1 = none.
    pub frz_pos: i64,
    /// frozen modified offset — `fa_frz_mod`. -1 = none.
    pub frz_mod: i64,
    pub prev_base: u8,
    pub prev_base_pos: i64,
    pub prev_is_insert: bool,
    /// TO_UPPER / TO_LOWER / -1 (unset)
    pub case: i8,
    pub napplied: u64,
    pub buf: Vec<u8>,
    /// region chromosome name (for mask overlap queries)
    pub chr_name: String,
}

impl ApplyState {
    pub fn new(ori_pos: i64, ref_seq: Vec<u8>) -> Self {
        ApplyState {
            ori_pos,
            mod_off: 0,
            frz_pos: -1,
            frz_mod: -1,
            prev_base: 0,
            prev_base_pos: -1,
            prev_is_insert: false,
            case: -1,
            napplied: 0,
            buf: ref_seq,
            chr_name: String::new(),
        }
    }

    pub fn with_chr(mut self, chr: impl Into<String>) -> Self {
        self.chr_name = chr.into();
        self
    }
}

/// Convenience: apply a sorted slice of records over a ref buffer, then run
/// the region-end absent fill and mask pass. `chr` is needed for the mask
/// overlap query. An optional `chain` accumulates ref↔alt gaps.
pub fn apply_region(
    chr: &str,
    ref_seq: Vec<u8>,
    ori_pos: i64,
    records: &[&VcfRecord],
    opts: &ApplyOptions,
    chain: Option<&mut Chain>,
) -> ApplyState {
    apply_region_planned(chr, ref_seq, ori_pos, records, opts, chain, None)
}

pub fn apply_region_planned(
    chr: &str,
    ref_seq: Vec<u8>,
    ori_pos: i64,
    records: &[&VcfRecord],
    opts: &ApplyOptions,
    chain: Option<&mut Chain>,
    plan: Option<&RegionPlan>,
) -> ApplyState {
    if records.is_empty() {
        return apply_empty_region(chr, ref_seq, ori_pos, opts);
    }

    if chain.is_none() {
        match plan.map(|p| p.lane) {
            Some(FastPathLane::SameLenOnly) => {
                if let Some(state) =
                    try_apply_same_len_region(chr, ref_seq.as_slice(), ori_pos, records, opts)
                {
                    return state;
                }
            }
            Some(FastPathLane::NormalizedEditScript | FastPathLane::MixedSimpleEdits) => {
                if let Some((state, _lane)) =
                    try_apply_edit_script_region(chr, ref_seq.as_slice(), ori_pos, records, opts)
                {
                    return state;
                }
            }
            Some(FastPathLane::FallbackStateMachine) => {}
            _ => {
                if let Some(state) =
                    try_apply_same_len_region(chr, ref_seq.as_slice(), ori_pos, records, opts)
                {
                    return state;
                }
                if let Some((state, _lane)) =
                    try_apply_edit_script_region(chr, ref_seq.as_slice(), ori_pos, records, opts)
                {
                    return state;
                }
            }
        }
    }

    let mut state = ApplyState::new(ori_pos, ref_seq).with_chr(chr);
    // bcftools masks each fa line on read, *before* any variant is applied
    // (consensus.c: mask_region runs in the fa-read loop; apply_variant runs
    // later on fa_buf). We mirror that here: mask the raw ref first, then apply
    // variants. Two consequences this fixes:
    //   1. Coordinate mapping: at this point buf is the raw ref, so original
    //      mask coords [mbeg,mend] map 1:1 to buf indices — no indel offset to
    //      account for (the old post-apply mask mislocated spans after indels).
    //   2. UC/LC case order: variants now see the masked case and sync to it,
    //      matching bcftools; the old post-apply mask overwrote variant case.
    if let Some(mask) = &opts.mask {
        mask.apply_to_buf(chr, &mut state.buf, state.ori_pos);
    }
    // hand the chain to apply_variant via a temporary Option<&mut> we can reborrow
    let mut chain_opt = chain;
    for rec in records {
        let ch = chain_opt.as_deref_mut();
        apply_variant(rec, &mut state, opts, ch);
    }
    // region-end absent fill (consensus.c:1177 apply_absent(HTS_POS_MAX))
    if let Some(absent) = opts.absent_allele {
        apply_absent(&mut state, i64::MAX, absent);
    }
    state
}

pub fn apply_region_with_stats(
    chr: &str,
    ref_seq: Vec<u8>,
    ori_pos: i64,
    records: &[&VcfRecord],
    opts: &ApplyOptions,
    stats: &mut RuntimeStats,
) -> ApplyState {
    stats.observe_region();
    for _ in records {
        stats.observe_record();
    }

    if records.is_empty() {
        stats.observe_lane(FastPathLane::EmptyRegion);
        return apply_empty_region(chr, ref_seq, ori_pos, opts);
    }

    if let Some(state) = try_apply_same_len_region(chr, ref_seq.as_slice(), ori_pos, records, opts)
    {
        stats.observe_lane(FastPathLane::SameLenOnly);
        for _ in 0..state.napplied {
            stats.observe_same_len_fastpath();
        }
        return state;
    }

    if let Some((state, lane)) =
        try_apply_edit_script_region(chr, ref_seq.as_slice(), ori_pos, records, opts)
    {
        stats.observe_lane(lane);
        for _ in 0..state.napplied {
            stats.observe_edit_script_fastpath();
        }
        return state;
    }

    let plan_opts = PlanOptions {
        chain: false,
        mark_del: opts.mark_del.is_some(),
        mark_ins: opts.mark_ins.is_some(),
        mark_snv: opts.mark_snv.is_some(),
        mask: opts.mask.is_some(),
        absent: opts.absent_allele.is_some(),
    };
    let plan = plan_region(records, plan_opts);
    stats.observe_lane(FastPathLane::FallbackStateMachine);
    stats.observe_fallback_records(records.len() as u64);
    if plan.fallback_reasons.is_empty() {
        stats.observe_fallback_reason(FallbackReason::UnsupportedMode);
    } else {
        for reason in plan.fallback_reasons {
            stats.observe_fallback_reason(reason);
        }
    }
    apply_region(chr, ref_seq, ori_pos, records, opts, None)
}

fn apply_empty_region(
    chr: &str,
    ref_seq: Vec<u8>,
    ori_pos: i64,
    opts: &ApplyOptions,
) -> ApplyState {
    let mut state = ApplyState::new(ori_pos, ref_seq).with_chr(chr);
    if let Some(mask) = &opts.mask {
        mask.apply_to_buf(chr, &mut state.buf, state.ori_pos);
    }
    if let Some(absent) = opts.absent_allele {
        apply_absent(&mut state, i64::MAX, absent);
    }
    state
}

struct SameLenPatch {
    idx: usize,
    pos: i64,
    rlen: usize,
    alt: AlleleBuf,
}

enum FastSelection {
    Skip,
    Allele {
        ialt: usize,
        alt_override: Option<AlleleBuf>,
    },
}

#[inline]
fn select_fastpath_allele(rec: &VcfRecord, opts: &ApplyOptions) -> Option<FastSelection> {
    if matches!(opts.sample_mode, SampleMode::ApplyAllAlt) {
        return Some(FastSelection::Allele {
            ialt: 1,
            alt_override: None,
        });
    }

    let selection = select_allele(rec, &opts.sample_mode, opts.missing_allele);
    match selection.ialt {
        None => Some(FastSelection::Skip),
        Some(i) if i >= 0 => Some(FastSelection::Allele {
            ialt: i as usize,
            alt_override: selection.alt_override,
        }),
        Some(_) => None,
    }
}

fn try_apply_same_len_region(
    chr: &str,
    ref_seq: &[u8],
    ori_pos: i64,
    records: &[&VcfRecord],
    opts: &ApplyOptions,
) -> Option<ApplyState> {
    if records.is_empty()
        || opts.absent_allele.is_some()
        || opts.missing_allele.is_some()
        || opts.mark_del.is_some()
        || opts.mark_ins.is_some()
        || opts.mark_snv.is_some()
        || opts.mask.is_some()
    {
        return None;
    }

    let mut patches: Vec<SameLenPatch> = Vec::with_capacity(records.len());
    let mut frz_pos = -1i64;
    for rec in records {
        if rec.alleles.len() == 1 {
            continue;
        }
        if rec.pos <= frz_pos {
            return None;
        }
        let (ialt, alt_override) = match select_fastpath_allele(rec, opts)? {
            FastSelection::Skip => continue,
            FastSelection::Allele { ialt, alt_override } => (ialt, alt_override),
        };
        if ialt >= rec.alleles.len() || rec.rlen <= 0 {
            return None;
        }
        let rlen = rec.rlen as usize;
        let idx = rec.pos - ori_pos;
        if idx < 0 {
            return None;
        }
        let idx = idx as usize;
        if idx + rlen > ref_seq.len() || rec.alleles[0].len() != rlen {
            return None;
        }
        if !ref_seq[idx..idx + rlen].eq_ignore_ascii_case(&rec.alleles[0]) {
            return None;
        }

        let alt = match alt_override {
            Some(alt) if alt.len() == rlen => alt,
            Some(_) => return None,
            None if rec.compiled.same_len_allele(ialt) => SmallVec::from_slice(&rec.alleles[ialt]),
            None => return None,
        };
        patches.push(SameLenPatch {
            idx,
            pos: rec.pos,
            rlen,
            alt,
        });
        frz_pos = rec.ref_end();
    }

    if patches.is_empty() {
        return None;
    }

    let mut state = ApplyState::new(ori_pos, ref_seq.to_vec()).with_chr(chr);
    for patch in patches {
        let first_base = state.buf[patch.idx];
        let last_base = state.buf[patch.idx + patch.rlen - 1];
        let to_upper = first_base.is_ascii_uppercase();
        state.case = if to_upper { TO_UPPER } else { TO_LOWER };
        copy_alt_with_case(
            &mut state.buf[patch.idx..patch.idx + patch.rlen],
            &patch.alt,
            to_upper,
        );
        state.prev_base = last_base;
        state.prev_base_pos = patch.pos + patch.rlen as i64 - 1;
        state.prev_is_insert = false;
        state.frz_mod = patch.idx as i64 + patch.rlen as i64;
        state.frz_pos = patch.pos + patch.rlen as i64 - 1;
        state.napplied += 1;
    }
    Some(state)
}

struct EditScriptPatch {
    idx: usize,
    pos: i64,
    rlen: usize,
    len_diff: i64,
    alt: AlleleBuf,
}

fn try_apply_edit_script_region(
    chr: &str,
    ref_seq: &[u8],
    ori_pos: i64,
    records: &[&VcfRecord],
    opts: &ApplyOptions,
) -> Option<(ApplyState, FastPathLane)> {
    if records.is_empty()
        || opts.absent_allele.is_some()
        || opts.missing_allele.is_some()
        || opts.mark_del.is_some()
        || opts.mark_ins.is_some()
        || opts.mark_snv.is_some()
        || opts.mask.is_some()
    {
        return None;
    }

    let mut patches: Vec<EditScriptPatch> = Vec::with_capacity(records.len());
    let mut frz_pos = -1i64;
    let mut total_delta = 0i64;
    let mut saw_same_len = false;
    let mut saw_len_change = false;

    for rec in records {
        if rec.alleles.len() == 1 {
            continue;
        }
        if rec.pos <= frz_pos {
            return None;
        }

        let (ialt, alt_override) = match select_fastpath_allele(rec, opts)? {
            FastSelection::Skip => continue,
            FastSelection::Allele { ialt, alt_override } => (ialt, alt_override),
        };
        if ialt >= rec.alleles.len() || rec.rlen <= 0 {
            return None;
        }
        let op = rec.compiled.allele_op(ialt)?;
        if !matches!(
            op.kind,
            AlleleOpKind::Ref | AlleleOpKind::SameLen | AlleleOpKind::Insert | AlleleOpKind::Delete
        ) {
            return None;
        }

        let rlen = rec.rlen as usize;
        let ref_allele = &rec.alleles[0];
        if ref_allele.len() != rlen || op.ref_len as usize != rlen {
            return None;
        }
        let idx = rec.pos - ori_pos;
        if idx < 0 {
            return None;
        }
        let idx = idx as usize;
        if idx + rlen > ref_seq.len() {
            return None;
        }
        if !ref_seq[idx..idx + rlen].eq_ignore_ascii_case(ref_allele) {
            return None;
        }

        let alt = match alt_override {
            Some(alt) if alt.len() == rlen && op.is_same_len_fastpath() => alt,
            Some(_) => return None,
            None => SmallVec::from_slice(&rec.alleles[ialt]),
        };
        if alt.starts_with(b"<") || alt.as_slice() == b"*" || alt.is_empty() {
            return None;
        }
        if op.alt_len as usize != alt.len() {
            return None;
        }
        if matches!(op.kind, AlleleOpKind::Insert | AlleleOpKind::Delete) && op.trim_beg != 1 {
            return None;
        }

        let len_diff = alt.len() as i64 - rlen as i64;
        if len_diff == 0 {
            saw_same_len = true;
        } else {
            saw_len_change = true;
        }
        total_delta += len_diff;
        patches.push(EditScriptPatch {
            idx,
            pos: rec.pos,
            rlen,
            len_diff,
            alt,
        });
        frz_pos = rec.ref_end();
    }

    if patches.is_empty() || !saw_len_change {
        return None;
    }
    let final_len = ref_seq.len() as i64 + total_delta;
    if final_len < 0 {
        return None;
    }

    let mut out = Vec::with_capacity(final_len as usize);
    let mut cursor = 0usize;
    let mut prev_base = 0u8;
    let mut prev_base_pos = -1i64;
    let mut prev_is_insert = false;
    let mut last_case = -1i8;
    let mut last_frz_pos = -1i64;
    let mut last_frz_mod = -1i64;

    for patch in &patches {
        if patch.idx < cursor || patch.idx + patch.rlen > ref_seq.len() {
            return None;
        }
        out.extend_from_slice(&ref_seq[cursor..patch.idx]);

        let first_base = ref_seq[patch.idx];
        let last_base = ref_seq[patch.idx + patch.rlen - 1];
        let to_upper = first_base.is_ascii_uppercase();
        last_case = if to_upper { TO_UPPER } else { TO_LOWER };
        extend_alt_with_case(&mut out, &patch.alt, to_upper);

        cursor = patch.idx + patch.rlen;
        prev_base = last_base;
        prev_base_pos = patch.pos + patch.rlen as i64 - 1;
        prev_is_insert = patch.len_diff > 0;
        last_frz_pos = prev_base_pos;
        last_frz_mod = out.len() as i64;
    }
    out.extend_from_slice(&ref_seq[cursor..]);

    let lane = if saw_same_len {
        FastPathLane::MixedSimpleEdits
    } else {
        FastPathLane::NormalizedEditScript
    };
    let mut state = ApplyState::new(ori_pos, out).with_chr(chr);
    state.mod_off = total_delta;
    state.frz_pos = last_frz_pos;
    state.frz_mod = last_frz_mod;
    state.prev_base = prev_base;
    state.prev_base_pos = prev_base_pos;
    state.prev_is_insert = prev_is_insert;
    state.case = last_case;
    state.napplied = patches.len() as u64;
    Some((state, lane))
}

// ---------------------------------------------------------------------------
// apply_absent (consensus.c:464) — fill gaps not covered by any record.
// ---------------------------------------------------------------------------

pub fn apply_absent(state: &mut ApplyState, pos: i64, absent: u8) {
    let blen = state.buf.len() as i64;
    if blen == 0 {
        return;
    }
    // if pos==frz+1, no gap (ie==ib); also no fill before region start.
    if pos <= state.frz_pos + 1 {
        return;
    }
    if pos <= state.ori_pos {
        return;
    }
    // bcftools: ie = (pos && pos-ori+off < blen) ? pos-ori+off : blen
    // Use saturating arithmetic: pos can be HTS_POS_MAX (i64::MAX) as a sentinel
    // for the region-end fill, in which case ie must clamp to blen.
    let cand = (pos - state.ori_pos).saturating_add(state.mod_off);
    let ie = if pos != 0 && cand < blen { cand } else { blen };
    let ib = if state.frz_mod < 0 { 0 } else { state.frz_mod };
    if ib < 0 || ie <= ib {
        return;
    }
    let ib = ib.max(0) as usize;
    let ie = (ie as usize).min(state.buf.len());
    state.buf[ib..ie].fill(absent);
}

// ---------------------------------------------------------------------------
// freeze_ref (consensus.c:474) — mark a ref block as frozen (no further apply).
// ---------------------------------------------------------------------------

fn freeze_ref(state: &mut ApplyState, rec: &VcfRecord) {
    if state.frz_pos >= rec.pos + rec.rlen as i64 - 1 {
        return;
    }
    state.frz_pos = rec.pos + rec.rlen as i64 - 1;
    state.frz_mod = rec.pos - state.ori_pos + state.mod_off + rec.rlen as i64;
}

// ---------------------------------------------------------------------------
// apply_variant (consensus.c:583)
// ---------------------------------------------------------------------------

/// Direct port of bcftools `apply_variant`. The `alt_len`/`rlen`/`alen`
/// assignments mirror the C control flow exactly; some intermediate writes are
/// overwritten on the next branch, which rustc flags — silence here.
#[allow(clippy::too_many_lines)]
#[allow(unused_assignments)]
pub fn apply_variant(
    rec: &VcfRecord,
    state: &mut ApplyState,
    opts: &ApplyOptions,
    chain: Option<&mut Chain>,
) {
    // 586: fill absent up to this record's pos (only with -a)
    if let Some(absent) = opts.absent_allele {
        apply_absent(state, rec.pos, absent);
    }
    // 760-765: ref-only record, nothing to apply (unless -M/-a which M2 doesn't set)
    if rec.alleles.len() == 1 && opts.missing_allele.is_none() && opts.absent_allele.is_none() {
        return;
    }
    // 590-600: mask check — char-mode mask overlapping this variant → skip.
    if let Some(mask) = &opts.mask {
        if mask.with.skips_variants() {
            // need chr for the query; bcftools uses args->chr (the region's chr).
            // We don't carry chr in state; the engine sets it via opts? For now
            // the mask check uses rec.rid-less name lookup through the engine.
            // M4: the engine passes the chr by storing it on ApplyState.
            if state.chr_name.as_bytes() != b""
                && mask.overlaps(&state.chr_name, rec.pos, rec.ref_end())
            {
                return;
            }
        }
    }

    // M3: allele selection via GT / -H / -I (consensus.c:602-758).
    let selection: AlleleSelection = select_allele(rec, &opts.sample_mode, opts.missing_allele);
    let ialt: i32 = match selection.ialt {
        None => return, // skip (overlap/missing-without-M/haplotype>ploidy warn)
        Some(i) => i,
    };

    // 766-788: missing allele (ialt == -1) — single pos or gvcf block.
    if ialt == -1 {
        if let Some(mchar) = opts.missing_allele {
            apply_missing(state, rec, mchar);
        }
        return;
    }

    // 760-765: ref-only record (n_allele==1) → freeze ref if -a, return.
    if rec.alleles.len() == 1 {
        if opts.absent_allele.is_some() {
            freeze_ref(state, rec);
        }
        return;
    }
    // Note: ialt==0 (REF selected) is NOT a return — bcftools applies REF
    // through the normal path (a no-op replacement that still updates frz_pos),
    // and the overlap trim below (consensus.c:811-825) may move it forward.
    let ialt_u = ialt as usize;
    if ialt_u >= rec.alleles.len() {
        // broken VCF (too few alts); bcftools errors. Skip defensively.
        return;
    }

    if chain.is_none()
        && try_apply_same_len_allele(rec, ialt_u, selection.alt_override.as_deref(), state, opts)
    {
        return;
    }

    // --- working copies of the mutable record fields (bcftools mutates rec) ---
    let mut pos = rec.pos;
    let mut rlen = rec.rlen as i64;

    // ref allele: front-trim offset into rec.alleles[0] (only increases).
    let ref_orig = &rec.alleles[0];
    let mut ref_off: i64 = 0;
    // alt allele: owned mutable copy. Use IUPAC-mixed bytes when select_allele
    // produced an override (the -I / iupac_GTs paths rewrite rec->d.allele[ialt]
    // in bcftools); otherwise copy the chosen allele verbatim.
    let mut alt_buf: AlleleBuf = match selection.alt_override {
        Some(v) => v,
        None => SmallVec::from_slice(&rec.alleles[ialt_u]),
    };
    let mut alt_off: i64 = 0;
    let mut alt_len: i64 = alt_buf.len() as i64;

    // 792-808: trim_beg / var_len / var_type
    let alt_is_symbolic_del = alt_buf.eq_ignore_ascii_case(b"<DEL>");
    let alt_is_symbolic_ins =
        alt_buf.len() >= 4 && alt_buf[0] == b'<' && alt_buf[1..4].eq_ignore_ascii_case(b"INS");
    let is_gvcf =
        alt_buf.eq_ignore_ascii_case(b"<*>") || alt_buf.eq_ignore_ascii_case(b"<NON_REF>");

    let is_indel = !alt_buf.starts_with(b"<") && alt_buf.as_slice() != b"*" && alt_len != rlen;
    let mut trim_beg: i64 = 0;
    let mut var_len: i64 = 0;
    if is_indel {
        // 794-797: first base same → anchor, trim it
        trim_beg = if ref_first(ref_orig, ref_off) == alt_first(&alt_buf, alt_off) {
            1
        } else {
            0
        };
        var_len = alt_len - rlen;
    } else if alt_is_symbolic_del {
        // 798-801
        trim_beg = 1;
        var_len = 1 - rlen;
    } else if alt_is_symbolic_ins {
        trim_beg = 1;
    }

    // 811-825: overlapping REF allele trim. When the REF allele was selected
    // (ialt==0) and it overlaps a previous deletion (pos <= frz_pos but ends
    // beyond it), trim the REF forward so it starts at frz_pos+1. bcftools
    // mutates rec->pos/rlen/allele[0] (REF and alt are the same string); we
    // mutate the local pos/rlen/ref_off, and advance alt_off too since alt_buf
    // is a copy of the REF allele in this path.
    if ialt == 0 && pos <= state.frz_pos && pos + rlen - 1 > state.frz_pos {
        let ntrim = state.frz_pos - pos + 1;
        let nref = ref_orig.len() as i64 - ref_off;
        if ntrim >= nref {
            // bcftools errors (unnormalized VCF); skip defensively.
            return;
        }
        pos += ntrim;
        rlen -= ntrim;
        ref_off += ntrim;
        alt_off += ntrim;
    }

    // 826-840: overlap skip check
    if pos <= state.frz_pos {
        // Can be OK iff this is an insertion not following another insertion (#888).
        let overlap = pos < state.frz_pos || trim_beg == 0 || var_len == 0 || state.prev_is_insert;
        if overlap {
            // bcftools prints a warning and skips; we skip silently (M2).
            return;
        }
    }

    // 847: idx = rec.pos - ori_pos + mod_off
    let mut idx: i64 = pos - state.ori_pos + state.mod_off;

    // 847-893: idx < 0 handling (variant starts before region but overlaps in)
    if idx < 0 {
        if alt_buf.starts_with(b"<") {
            // symbolic: shift so anchor sits at idx=-1
            pos -= idx + 1;
            rlen += idx + 1;
            idx = -1;
        } else if (ref_orig.len() as i64 - ref_off) < -idx {
            // ref allele shorter than the overhang — shouldn't happen; skip
            return;
        } else if alt_len - alt_off > -idx {
            // ref and alt both overlap the fa buffer
            pos -= idx;
            rlen += idx;
            ref_off -= idx; // -= idx (idx<0) → advances forward
            alt_off -= idx;
            alt_len = alt_buf.len() as i64 - alt_off; // recompute strlen
            idx = 0;
        } else {
            // ref overlaps fa but alt does not: trim to leave one base before
            pos -= idx + 1;
            rlen += idx + 1;
            ref_off -= idx + 1;
            // alt_allele += strlen(alt_allele)-1  → point to last char
            alt_off = alt_off + alt_len - 1;
            alt_len = 1;
            idx = -1;
        }
    }

    // 870-882: rlen exceeds available buffer → trim
    let blen = state.buf.len() as i64;
    if rlen > blen - idx {
        rlen = blen - idx;
        if !alt_buf.starts_with(b"<") {
            let alen = alt_len - alt_off;
            if alen > rlen {
                alt_len = alt_off + rlen; // truncate alt to rlen
            }
        }
    }

    // 884: variant entirely beyond the buffer
    if idx > 0 && idx >= blen {
        return;
    }

    // 896-924: symbolic allele handling
    let mut len_diff: i64;
    let mut alen: i64;
    if alt_buf.starts_with(b"<") {
        // bcftools (consensus.c:899): only <DEL>, <*>, <NON_REF> are supported.
        // Any other symbolic allele (<INS>, <DUP>, ...) → bcftools errors; we
        // skip defensively. Note <INS> sets trim_beg above (line 808) but is
        // still rejected here, matching bcftools.
        if !alt_is_symbolic_del && !is_gvcf {
            return;
        }
        if alt_is_symbolic_del {
            if opts.mark_del.is_some() {
                // M4: mark_del path (placeholder — M4 implements fully)
                alt_buf = mark_del_bytes(ref_orig, ref_off, rlen, None, opts.mark_del);
                alt_off = 0;
                alt_len = rlen;
                alen = rlen;
                len_diff = 0;
            } else {
                len_diff = 1 - rlen;
                // alt_allele = ref_allele; alen = 1
                alt_buf = SmallVec::from_slice(&ref_orig[ref_off as usize..]);
                alt_off = 0;
                alt_len = 1;
                alen = 1;
            }
        } else {
            // <*> or <NON_REF> — gVCF reference block: freeze and skip
            freeze_ref(state, rec);
            return;
        }
    } else if idx >= 0 && !ref_matches(state, idx, ref_orig, ref_off, rlen) {
        // 925-965: REF mismatch with the fa buffer — prev_base fallback or error
        let fail = ref_mismatch_fallback_ok(state, pos, rlen, ref_orig, ref_off, idx);
        if !fail {
            // bcftools errors here. We surface an error by skipping (M2).
            return;
        }
        alen = alt_len - alt_off;
        len_diff = alen - rlen;
        if opts.mark_del.is_some() && len_diff < 0 {
            alt_buf = mark_del_bytes(
                ref_orig,
                ref_off,
                rlen,
                Some(&alt_buf[alt_off as usize..alt_off as usize + alt_len as usize]),
                opts.mark_del,
            );
            alt_off = 0;
            alt_len = rlen;
            alen = rlen;
            len_diff = 0;
        }
    } else {
        // 967-976: ref matches
        alen = alt_len - alt_off;
        len_diff = alen - rlen;
        if opts.mark_del.is_some() && len_diff < 0 {
            alt_buf = mark_del_bytes(
                ref_orig,
                ref_off,
                rlen,
                Some(&alt_buf[alt_off as usize..alt_off as usize + alt_len as usize]),
                opts.mark_del,
            );
            alt_off = 0;
            alt_len = rlen;
            alen = rlen;
            len_diff = 0;
        }
    }

    // 992-996: case sync
    let safe_idx = if idx < 0 { 0 } else { idx as usize };
    let base = state.buf[safe_idx.min(state.buf.len() - 1)];
    state.case = if base.is_ascii_uppercase() {
        TO_UPPER
    } else {
        TO_LOWER
    };
    let apply_case = |b: u8| -> u8 {
        if state.case == TO_UPPER {
            b.to_ascii_uppercase()
        } else {
            b.to_ascii_lowercase()
        }
    };
    for i in 0..alen as usize {
        let p = (alt_off as usize) + i;
        if p < alt_buf.len() {
            alt_buf[p] = apply_case(alt_buf[p]);
        }
    }

    // 998-1001: mark_ins (len_diff>0) / mark_snv.
    if let Some(mark) = opts.mark_ins {
        if len_diff > 0 {
            mark_ins_bytes(ref_orig, ref_off, &mut alt_buf, alt_off, alen, mark);
        }
    }
    if let Some(mark) = opts.mark_snv {
        mark_snv_bytes(ref_orig, ref_off, &mut alt_buf, alt_off, alen, mark);
    }

    // 1003-1040: the actual replacement
    if len_diff <= 0 {
        // deletion or same-size
        // prev_base = buf[idx+rlen-1]
        let pb_idx = idx + rlen - 1;
        if pb_idx >= 0 && (pb_idx as usize) < state.buf.len() {
            state.prev_base = state.buf[pb_idx as usize];
        }
        state.prev_base_pos = pos + rlen - 1;
        state.prev_is_insert = false;
        state.frz_mod = idx + alen;

        // write alt[trim_beg..alen] into buf[idx..]
        for i in trim_beg..alen {
            let wp = idx + i;
            if wp >= 0 && (wp as usize) < state.buf.len() {
                let ap = (alt_off as usize) + i as usize;
                if ap < alt_buf.len() {
                    state.buf[wp as usize] = alt_buf[ap];
                }
            }
        }
        // memmove shrink: remove buf[idx+alen .. idx+rlen]
        let drain_start = (idx + alen).max(0) as usize;
        let drain_end = (idx + rlen).max(0) as usize;
        if drain_end > drain_start && drain_end <= state.buf.len() {
            remove_range(&mut state.buf, drain_start, drain_end);
        }
    } else {
        // insertion (len_diff > 0)
        state.prev_is_insert = true;
        state.prev_base_pos = pos;

        // make room: insert len_diff bytes at buf[idx+rlen]
        let ins_pos = ((idx + rlen).max(0) as usize).min(state.buf.len());
        insert_zeros(&mut state.buf, ins_pos, len_diff as usize);

        // ibeg: skip leading matching bases (anchor) unchanged by the insertion,
        // so we don't overwrite a preceding variant's base (#888).
        let mut ibeg: i64 = 0;
        while ibeg < alen {
            let rb = ref_off as usize + ibeg as usize;
            let ab = alt_off as usize + ibeg as usize;
            if rb >= ref_orig.len() || ab >= alt_buf.len() {
                break;
            }
            if ref_orig[rb] != alt_buf[ab] {
                break;
            }
            if pos + ibeg > state.prev_base_pos {
                break;
            }
            ibeg += 1;
        }
        for i in ibeg..alen {
            let wp = idx + i;
            if wp >= 0 && (wp as usize) < state.buf.len() {
                let ap = (alt_off as usize) + i as usize;
                if ap < alt_buf.len() {
                    state.buf[wp as usize] = alt_buf[ap];
                }
            }
        }
        state.frz_mod = idx + alen - ibeg + 1;
    }

    // 1041-1054: chain
    if let Some(chain) = chain {
        if len_diff != 0 {
            // ref/alt allele first-base views for the anchor check
            let ref_b0 = ref_first(ref_orig, ref_off);
            let alt_b0 = alt_first(&alt_buf, alt_off);
            if ascii_eq_ignore_case(ref_b0, alt_b0) {
                // first base same → extend block by 1bp: start+1, alleles 1bp shorter
                chain.push_gap(pos + 1, rlen - 1, pos + 1 + state.mod_off, alen - 1);
            } else {
                chain.push_gap(pos, rlen, pos + state.mod_off, alen);
            }
        }
    }

    // 1055-1058: state update
    // buf.l += len_diff  (Vec len already updated via drain/splice)
    state.mod_off += len_diff;
    state.frz_pos = pos + rlen - 1;
    state.napplied += 1;
}

// ---------------------------------------------------------------------------
// small helpers
// ---------------------------------------------------------------------------

/// Missing-GT path (consensus.c:766-788): replace the REF span with the missing
/// char. For a gvcf reference block (rlen>1, var_type is REF) the whole block
/// is overwritten; otherwise a single position. bcftools achieves this by
/// rewriting alleles to `REF,missing_char` and setting ialt=1; we inline the
/// same buffer effect.
fn apply_missing(state: &mut ApplyState, rec: &VcfRecord, mchar: u8) {
    let idx = rec.pos - state.ori_pos + state.mod_off;
    if idx < 0 || idx as usize >= state.buf.len() {
        return;
    }
    let rlen = rec.rlen as usize;
    let end = ((idx as usize) + rlen).min(state.buf.len());
    state.buf[idx as usize..end].fill(mchar);
    // freeze so subsequent absent fill / overlap checks see this position done
    state.frz_pos = rec.pos + rec.rlen as i64 - 1;
    state.frz_mod = idx + (end - idx as usize) as i64;
    state.prev_base = mchar;
    state.prev_base_pos = rec.pos + rec.rlen as i64 - 1;
    state.prev_is_insert = false;
    state.napplied += 1;
}

fn ref_first(ref_allele: &[u8], ref_off: i64) -> u8 {
    let i = ref_off as usize;
    if i < ref_allele.len() {
        ref_allele[i]
    } else {
        0
    }
}

fn alt_first(alt_buf: &[u8], alt_off: i64) -> u8 {
    let i = alt_off as usize;
    if i < alt_buf.len() {
        alt_buf[i]
    } else {
        0
    }
}

fn try_apply_same_len_allele(
    rec: &VcfRecord,
    ialt: usize,
    alt_override: Option<&[u8]>,
    state: &mut ApplyState,
    opts: &ApplyOptions,
) -> bool {
    if opts.mark_del.is_some() || opts.mark_ins.is_some() || opts.mark_snv.is_some() {
        return false;
    }
    if rec.rlen <= 0 {
        return false;
    }
    let rlen = rec.rlen as usize;
    let ref_allele = &rec.alleles[0];
    let alt_allele = match alt_override {
        Some(a) => {
            if a.len() != rlen {
                return false;
            }
            a
        }
        None => match rec.alleles.get(ialt) {
            Some(a) if rec.compiled.same_len_allele(ialt) => &a[..],
            Some(_) => return false,
            None => return false,
        },
    };
    if ref_allele.len() != rlen
        || alt_allele.len() != rlen
        || alt_allele.starts_with(b"<")
        || alt_allele == b"*"
    {
        return false;
    }

    let pos = rec.pos;
    if pos <= state.frz_pos {
        if ialt == 0 && pos + rec.rlen as i64 - 1 > state.frz_pos {
            return false;
        }
        return true;
    }

    let idx = pos - state.ori_pos + state.mod_off;
    if idx < 0 {
        return false;
    }
    let idx = idx as usize;
    if idx >= state.buf.len() {
        return true;
    }
    if idx + rlen > state.buf.len() {
        return false;
    }

    if !state.buf[idx..idx + rlen].eq_ignore_ascii_case(ref_allele)
        && !ref_mismatch_fallback_ok(state, pos, rec.rlen as i64, ref_allele, 0, idx as i64)
    {
        return true;
    }

    let first_base = state.buf[idx];
    let last_base = state.buf[idx + rlen - 1];
    let to_upper = first_base.is_ascii_uppercase();
    state.case = if to_upper { TO_UPPER } else { TO_LOWER };
    copy_alt_with_case(&mut state.buf[idx..idx + rlen], alt_allele, to_upper);

    state.prev_base = last_base;
    state.prev_base_pos = pos + rec.rlen as i64 - 1;
    state.prev_is_insert = false;
    state.frz_mod = idx as i64 + rec.rlen as i64;
    state.frz_pos = pos + rec.rlen as i64 - 1;
    state.napplied += 1;
    true
}

#[inline]
fn ascii_eq_ignore_case(a: u8, b: u8) -> bool {
    a.eq_ignore_ascii_case(&b)
}

#[inline]
fn needs_uppercase_conversion(bytes: &[u8]) -> bool {
    bytes.iter().any(u8::is_ascii_lowercase)
}

#[inline]
fn needs_lowercase_conversion(bytes: &[u8]) -> bool {
    bytes.iter().any(u8::is_ascii_uppercase)
}

fn copy_alt_with_case(dst: &mut [u8], alt: &[u8], to_upper: bool) {
    debug_assert_eq!(dst.len(), alt.len());
    if alt.len() == 1 {
        dst[0] = if to_upper {
            alt[0].to_ascii_uppercase()
        } else {
            alt[0].to_ascii_lowercase()
        };
        return;
    }
    let needs_conversion = if to_upper {
        needs_uppercase_conversion(alt)
    } else {
        needs_lowercase_conversion(alt)
    };
    if !needs_conversion {
        dst.copy_from_slice(alt);
        return;
    }
    if to_upper {
        for (d, &src) in dst.iter_mut().zip(alt) {
            *d = src.to_ascii_uppercase();
        }
    } else {
        for (d, &src) in dst.iter_mut().zip(alt) {
            *d = src.to_ascii_lowercase();
        }
    }
}

fn extend_alt_with_case(out: &mut Vec<u8>, alt: &[u8], to_upper: bool) {
    if alt.len() == 1 {
        out.push(if to_upper {
            alt[0].to_ascii_uppercase()
        } else {
            alt[0].to_ascii_lowercase()
        });
        return;
    }
    let needs_conversion = if to_upper {
        needs_uppercase_conversion(alt)
    } else {
        needs_lowercase_conversion(alt)
    };
    if !needs_conversion {
        out.extend_from_slice(alt);
        return;
    }
    if to_upper {
        out.extend(alt.iter().map(u8::to_ascii_uppercase));
    } else {
        out.extend(alt.iter().map(u8::to_ascii_lowercase));
    }
}

fn remove_range(buf: &mut Vec<u8>, start: usize, end: usize) {
    if end <= start || start >= buf.len() {
        return;
    }
    let end = end.min(buf.len());
    let len = buf.len();
    if end < len {
        buf.copy_within(end..len, start);
    }
    buf.truncate(len - (end - start));
}

fn insert_zeros(buf: &mut Vec<u8>, pos: usize, n: usize) {
    if n == 0 {
        return;
    }
    let pos = pos.min(buf.len());
    let old_len = buf.len();
    buf.resize(old_len + n, 0);
    if pos < old_len {
        buf.copy_within(pos..old_len, pos + n);
        buf[pos..pos + n].fill(0);
    }
}

/// `mark_ins` (consensus.c:502): for the inserted bases (alt[nref..nalt]) apply
/// the mark char / case. `alen` is the effective alt length, `alt_off` its start.
fn mark_ins_bytes(
    ref_allele: &[u8],
    ref_off: i64,
    alt_buf: &mut [u8],
    alt_off: i64,
    alen: i64,
    mark: u8,
) {
    let nref = (ref_allele.len() as i64 - ref_off).max(0);
    let start = (alt_off + nref) as usize;
    let end = (alt_off + alen) as usize;
    if start >= end || start >= alt_buf.len() {
        return;
    }
    let end = end.min(alt_buf.len());
    if mark == TO_UPPER as u8 {
        for b in &mut alt_buf[start..end] {
            *b = b.to_ascii_uppercase();
        }
    } else if mark == TO_LOWER as u8 {
        for b in &mut alt_buf[start..end] {
            *b = b.to_ascii_lowercase();
        }
    } else {
        for b in &mut alt_buf[start..end] {
            *b = mark;
        }
    }
}

/// `mark_snv` (consensus.c:511): for the overlapping ref/alt prefix where bases
/// differ, apply the mark char / case.
fn mark_snv_bytes(
    ref_allele: &[u8],
    ref_off: i64,
    alt_buf: &mut [u8],
    alt_off: i64,
    alen: i64,
    mark: u8,
) {
    let nref = (ref_allele.len() as i64 - ref_off).max(0);
    let n = nref.min(alen) as usize;
    let rstart = ref_off as usize;
    let astart = alt_off as usize;
    if rstart + n > ref_allele.len() || astart + n > alt_buf.len() {
        return;
    }
    if mark == TO_UPPER as u8 {
        for i in 0..n {
            if !ascii_eq_ignore_case(ref_allele[rstart + i], alt_buf[astart + i]) {
                alt_buf[astart + i] = alt_buf[astart + i].to_ascii_uppercase();
            }
        }
    } else if mark == TO_LOWER as u8 {
        for i in 0..n {
            if !ascii_eq_ignore_case(ref_allele[rstart + i], alt_buf[astart + i]) {
                alt_buf[astart + i] = alt_buf[astart + i].to_ascii_lowercase();
            }
        }
    } else {
        for i in 0..n {
            if !ascii_eq_ignore_case(ref_allele[rstart + i], alt_buf[astart + i]) {
                alt_buf[astart + i] = mark;
            }
        }
    }
}

/// `strncasecmp(ref_allele, fa_buf+idx, rlen)` (consensus.c:925)
fn ref_matches(state: &ApplyState, idx: i64, ref_allele: &[u8], ref_off: i64, rlen: i64) -> bool {
    let start = idx as usize;
    if start + rlen as usize > state.buf.len() {
        return false;
    }
    let rstart = ref_off as usize;
    if rstart + rlen as usize > ref_allele.len() {
        return false;
    }
    let n = rlen as usize;
    state.buf[start..start + n].eq_ignore_ascii_case(&ref_allele[rstart..rstart + n])
}

/// prev_base fallback (consensus.c:942-947): ok iff prev_base_pos==pos and
/// the ref's first base matches prev_base, with either rlen==1 or the rest
/// matching buf[idx+1..].
///
/// Uses the *local* `pos`/`rlen` (after idx<0 handling and overlap trim may
/// have modified them), matching bcftools which mutates rec->pos/rec->rlen
/// in place before reaching this check.
fn ref_mismatch_fallback_ok(
    state: &ApplyState,
    pos: i64,
    rlen: i64,
    ref_allele: &[u8],
    ref_off: i64,
    idx: i64,
) -> bool {
    if state.prev_base_pos != pos {
        return false;
    }
    let r0 = ref_first(ref_allele, ref_off);
    if !ascii_eq_ignore_case(r0, state.prev_base) {
        return false;
    }
    if rlen == 1 {
        return true;
    }
    // ref[1..rlen] == buf[idx+1..idx+rlen]
    let start = (idx + 1) as usize;
    let rstart = (ref_off as usize) + 1;
    let n = (rlen as usize) - 1;
    if start + n > state.buf.len() || rstart + n > ref_allele.len() {
        return false;
    }
    state.buf[start..start + n].eq_ignore_ascii_case(&ref_allele[rstart..rstart + n])
}

/// mark_del (consensus.c:482): produce an alt of length `rlen` where the first
/// nalt bytes come from `alt` (or ref for symbolic <DEL>) and the rest is `mark`.
fn mark_del_bytes(
    ref_allele: &[u8],
    ref_off: i64,
    rlen: i64,
    alt: Option<&[u8]>,
    mark: Option<u8>,
) -> AlleleBuf {
    let mark = match mark {
        Some(m) => m,
        None => return AlleleBuf::new(),
    };
    let mut out = AlleleBuf::with_capacity(rlen as usize);
    if let Some(alt) = alt {
        out.extend_from_slice(alt);
    } else {
        // symbolic <DEL>: copy ref
        let rstart = ref_off as usize;
        let nref = (ref_allele.len() - rstart).min(rlen as usize);
        out.extend_from_slice(&ref_allele[rstart..rstart + nref]);
    }
    while (out.len() as i64) < rlen {
        out.push(mark);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf_store::VcfStore;

    fn write_vcf(name: &str, body: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("consensus_rs_apply_{}", name));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let vcf = dir.join("t.vcf");
        let header = "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000>\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        std::fs::write(&vcf, format!("{}{}", header, body)).unwrap();
        vcf
    }

    fn apply_all(vcf_path: &std::path::Path, ref_seq: &[u8], ori_pos: i64, end: i64) -> Vec<u8> {
        let store = VcfStore::load(vcf_path).unwrap();
        let recs = store.query("chr1", ori_pos, end, 1);
        let opts = ApplyOptions::default();
        let mut state = ApplyState::new(ori_pos, ref_seq.to_vec());
        for r in recs {
            apply_variant(r, &mut state, &opts, None);
        }
        state.buf
    }

    /// ref = "ACGTACGT" (8 bp), 0-based positions 0..7
    const REF: &[u8] = b"ACGTACGT";

    #[test]
    fn snp_replaces_single_base() {
        // chr1:2 (1-based) = 0-based pos 1, REF=C ALT=G -> replace buf[1] C->G
        let vcf = write_vcf("snp", "chr1\t2\t.\tC\tG\t.\t.\t.\n");
        let out = apply_all(&vcf, REF, 0, 7);
        assert_eq!(out, b"AGGTACGT");
    }

    #[test]
    fn mnp_replaces_block() {
        // 0-based pos 2, REF=GT ALT=CA -> replace buf[2..4] "GT"->"CA"
        let vcf = write_vcf("mnp", "chr1\t3\t.\tGT\tCA\t.\t.\t.\n");
        let out = apply_all(&vcf, REF, 0, 7);
        assert_eq!(out, b"ACCAACGT");
    }

    #[test]
    fn stats_entry_uses_same_len_region_fastpath() {
        let vcf = write_vcf(
            "same_len_stats",
            "chr1\t2\t.\tC\tG\t.\t.\t.\n\
             chr1\t5\t.\tA\tT\t.\t.\t.\n",
        );
        let store = VcfStore::load(&vcf).unwrap();
        let recs = store.query("chr1", 0, 7, 1);
        let opts = ApplyOptions::default();
        let mut stats = RuntimeStats::default();

        let state = apply_region_with_stats("chr1", REF.to_vec(), 0, &recs, &opts, &mut stats);

        assert_eq!(state.buf, b"AGGTTCGT");
        assert_eq!(stats.regions_total, 1);
        assert_eq!(stats.records_seen, 2);
        assert_eq!(stats.lane_count(FastPathLane::SameLenOnly), 1);
        assert_eq!(stats.same_len_fastpath_records, 2);
        assert_eq!(stats.fallback_records, 0);
    }

    #[test]
    fn same_len_fastpath_syncs_lowercase_ref_case() {
        let vcf = write_vcf("same_len_lower", "chr1\t2\t.\tC\tG\t.\t.\t.\n");
        let store = VcfStore::load(&vcf).unwrap();
        let recs = store.query("chr1", 0, 7, 1);
        let opts = ApplyOptions::default();
        let mut stats = RuntimeStats::default();

        let state =
            apply_region_with_stats("chr1", b"acgtacgt".to_vec(), 0, &recs, &opts, &mut stats);

        assert_eq!(state.buf, b"aggtacgt");
        assert_eq!(stats.lane_count(FastPathLane::SameLenOnly), 1);
    }

    #[test]
    fn stats_entry_uses_empty_region_lane() {
        let opts = ApplyOptions::default();
        let mut stats = RuntimeStats::default();
        let empty: Vec<&VcfRecord> = Vec::new();

        let state = apply_region_with_stats("chr1", REF.to_vec(), 0, &empty, &opts, &mut stats);

        assert_eq!(state.buf, REF);
        assert_eq!(stats.regions_total, 1);
        assert_eq!(stats.records_seen, 0);
        assert_eq!(stats.lane_count(FastPathLane::EmptyRegion), 1);
        assert_eq!(stats.fallback_records, 0);
    }

    #[test]
    fn stats_entry_uses_normalized_edit_script_fastpath() {
        let vcf = write_vcf(
            "edit_script_stats",
            "chr1\t2\t.\tC\tCAA\t.\t.\t.\n\
             chr1\t5\t.\tACG\tA\t.\t.\t.\n",
        );
        let store = VcfStore::load(&vcf).unwrap();
        let recs = store.query("chr1", 0, 7, 1);
        let opts = ApplyOptions::default();
        let mut stats = RuntimeStats::default();

        let state = apply_region_with_stats("chr1", REF.to_vec(), 0, &recs, &opts, &mut stats);

        assert_eq!(state.buf, b"ACAAGTAT");
        assert_eq!(stats.lane_count(FastPathLane::NormalizedEditScript), 1);
        assert_eq!(stats.edit_script_fastpath_records, 2);
        assert_eq!(stats.fallback_records, 0);
    }

    #[test]
    fn edit_script_fastpath_syncs_lowercase_ref_case() {
        let vcf = write_vcf(
            "edit_script_lower",
            "chr1\t2\t.\tC\tCAA\t.\t.\t.\n\
             chr1\t5\t.\tACG\tA\t.\t.\t.\n",
        );
        let store = VcfStore::load(&vcf).unwrap();
        let recs = store.query("chr1", 0, 7, 1);
        let opts = ApplyOptions::default();
        let mut stats = RuntimeStats::default();

        let state =
            apply_region_with_stats("chr1", b"acgtacgt".to_vec(), 0, &recs, &opts, &mut stats);

        assert_eq!(state.buf, b"acaagtat");
        assert_eq!(stats.lane_count(FastPathLane::NormalizedEditScript), 1);
    }

    #[test]
    fn pure_insertion_grows_buffer() {
        // 0-based pos 3, REF=T ALT=TGGG -> insert "GGG" after T (anchor T kept)
        let vcf = write_vcf("ins", "chr1\t4\t.\tT\tTGGG\t.\t.\t.\n");
        let out = apply_all(&vcf, REF, 0, 7);
        assert_eq!(out, b"ACGTGGGACGT");
    }

    #[test]
    fn pure_deletion_shrinks_buffer() {
        // 0-based pos 3, REF=TAC ALT=T -> delete "AC" after anchor T
        let vcf = write_vcf("del", "chr1\t4\t.\tTAC\tT\t.\t.\t.\n");
        let out = apply_all(&vcf, REF, 0, 7);
        assert_eq!(out, b"ACGTGT");
    }

    #[test]
    fn complex_indel_insert_plus_delete() {
        // 0-based pos 3, REF=TAC ALT=TG -> anchor T, then "AC"->"G" (del 1, ins 1 net)
        let vcf = write_vcf("complex", "chr1\t4\t.\tTAC\tTG\t.\t.\t.\n");
        let out = apply_all(&vcf, REF, 0, 7);
        assert_eq!(out, b"ACGTGGT");
    }

    #[test]
    fn symbolic_del() {
        // <DEL> at 0-based pos 3, REF=TAC rlen=3 -> keep anchor T, delete "AC"
        let vcf = write_vcf("symboldel", "chr1\t4\t.\tTAC\t<DEL>\t.\t.\t.\n");
        let out = apply_all(&vcf, REF, 0, 7);
        assert_eq!(out, b"ACGTGT");
    }

    /// M2 acceptance: byte-for-byte parity with `bcftools consensus` (no sample,
    /// no -I). Ignored by default (needs bgzip/tabix/bcftools on PATH); run with
    /// `cargo test -- --ignored bcftools_parity`.
    #[test]
    #[ignore]
    fn bcftools_parity() {
        use std::process::Command;
        let dir = std::env::temp_dir().join("consensus_rs_apply_parity");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        // 60-bp ref so wrap behaviour is also exercised. Cycle ACGT.
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let seqb = seq.as_bytes();
        std::fs::write(dir.join("ref.fa"), format!(">chr1\n{}\n", seq)).unwrap();
        let ref_fa = dir.join("ref.fa");

        // 1-based REF base at position p1
        let r = |p1: usize| -> char { seqb[p1 - 1] as char };

        // Build a multi-variant VCF whose REF is read from the actual sequence,
        // so we never get coordinate/REF mismatches: SNP, insertion, deletion.
        let multi = format!(
            "chr1\t2\t.\t{r2}\tG\t.\t.\t.\n\
             chr1\t10\t.\t{r10}\t{r10}TT\t.\t.\t.\n\
             chr1\t12\t.\t{r12}{r13}{r14}\t{r12}\t.\t.\t.\n",
            r2 = r(2),
            r10 = r(10),
            r12 = r(12),
            r13 = r(13),
            r14 = r(14),
        );

        // insert-follows-deletion: v1 deletes "AC" (1-based 5-6), v2 inserts
        // "AA" after the base that now follows. v2's REF no longer matches the
        // modified buf → exercises the prev_base fallback (consensus.c:942-947).
        let ins_after_del = format!(
            "chr1\t4\t.\t{r4}{r5}{r6}\t{r4}\t.\t.\t.\n\
             chr1\t6\t.\t{r6}{r7}\t{r6}{r7}AA\t.\t.\t.\n",
            r4 = r(4),
            r5 = r(5),
            r6 = r(6),
            r7 = r(7),
        );

        let cases: Vec<(&str, &str)> = vec![
            ("snp", "chr1\t2\t.\tC\tG\t.\t.\t.\n"),
            ("mnp", "chr1\t3\t.\tGT\tCA\t.\t.\t.\n"),
            ("ins", "chr1\t4\t.\tT\tTGGG\t.\t.\t.\n"),
            ("del", "chr1\t4\t.\tTAC\tT\t.\t.\t.\n"),
            ("complex", "chr1\t4\t.\tTAC\tTG\t.\t.\t.\n"),
            ("symboldel", "chr1\t4\t.\tTAC\t<DEL>\t.\t.\t.\n"),
            ("multi", multi.as_str()),
            ("ins_after_del", ins_after_del.as_str()),
        ];

        for (name, body) in &cases {
            // write + bgzip + tabix
            let vcf = dir.join(format!("{}.vcf", name));
            let header = "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100>\n\
                #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
            std::fs::write(&vcf, format!("{}{}", header, body)).unwrap();
            let st = Command::new("bgzip").arg(&vcf).status().expect("bgzip");
            assert!(st.success());
            let gz = vcf.with_extension("vcf.gz");
            let st = Command::new("tabix")
                .arg("-p")
                .arg("vcf")
                .arg(&gz)
                .status()
                .expect("tabix");
            assert!(st.success());

            // ground truth: bcftools consensus (no -s, no -I)
            let out = Command::new("bcftools")
                .args(["consensus", "-f"])
                .arg(&ref_fa)
                .arg(&gz)
                .output()
                .expect("bcftools");
            assert!(
                out.status.success(),
                "bcftools failed for `{}`: {}",
                name,
                String::from_utf8_lossy(&out.stderr)
            );
            let theirs: Vec<u8> = String::from_utf8_lossy(&out.stdout)
                .lines()
                .filter(|l| !l.starts_with('>') && !l.is_empty())
                .flat_map(|l| l.as_bytes().iter().copied())
                .collect();

            // ours: apply over the whole ref region
            let store = VcfStore::load(&gz).unwrap();
            let end = seq.len() as i64 - 1;
            let recs = store.query("chr1", 0, end, 1);
            let opts = ApplyOptions::default();
            let mut state = ApplyState::new(0, seq.as_bytes().to_vec());
            for r in recs {
                apply_variant(r, &mut state, &opts, None);
            }
            assert_eq!(
                state.buf,
                theirs,
                "case `{}` mismatch:\n  ours:   {}\n  theirs: {}",
                name,
                String::from_utf8_lossy(&state.buf),
                String::from_utf8_lossy(&theirs)
            );
        }
    }

    /// M3 acceptance: parity with `bcftools consensus -s X -H ...` and `-I`.
    /// Ignored by default (needs bgzip/tabix/bcftools). Run with
    /// `cargo test -- --ignored bcftools_haplotype_parity`.
    #[test]
    #[ignore]
    fn bcftools_haplotype_parity() {
        use crate::haplotype::{HaplotypeSpec, SampleMode};
        use std::process::Command;
        let dir = std::env::temp_dir().join("consensus_rs_hap_parity");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let seqb = seq.as_bytes();
        std::fs::write(dir.join("ref.fa"), format!(">chr1\n{}\n", seq)).unwrap();
        let ref_fa = dir.join("ref.fa");
        let r = |p1: usize| -> char { seqb[p1 - 1] as char };

        // Two samples, diploid, mixed phased/unphased + a missing GT.
        // pos 2: C/G  het (S1=0|1 phased, S2=1/1)
        // pos 10: C>TT ins (S1=0/1, S2=./. missing)
        // pos 20: T>A snp (S1=1|0, S2=0/1 unphased -> IUPAC)
        let body = format!(
            "##fileformat=VCFv4.3\n\
             ##contig=<ID=chr1,length=100>\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"G\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
             chr1\t2\t.\t{r2}\tG\t.\t.\t.\tGT\t0|1\t1/1\n\
             chr1\t10\t.\t{r10}\t{r10}TT\t.\t.\t.\tGT\t0/1\t./.\n\
             chr1\t20\t.\t{r20}\tA\t.\t.\t.\tGT\t1|0\t0/1\n",
            r2 = r(2),
            r10 = r(10),
            r20 = r(20),
        );
        let vcf = dir.join("hap.vcf");
        std::fs::write(&vcf, body).unwrap();
        let st = Command::new("bgzip").arg(&vcf).status().expect("bgzip");
        assert!(st.success());
        let gz = vcf.with_extension("vcf.gz");
        let st = Command::new("tabix")
            .arg("-p")
            .arg("vcf")
            .arg(&gz)
            .status()
            .expect("tabix");
        assert!(st.success());

        let store = VcfStore::load(&gz).unwrap();
        let end = seq.len() as i64 - 1;
        let recs = store.query("chr1", 0, end, 1);

        // (cli_args, sample_mode). sample S1 = idx 0, S2 = idx 1.
        let modes: Vec<(Vec<&str>, SampleMode)> = vec![
            // -I (no sample): IUPAC mix REF+ALT
            (vec!["-I"], SampleMode::IupacFromRefAlt),
            // -s S1 -H R / A / I / 1 / 2 / 1pIu / 2pIu / LR / LA / SR / SA
            (
                vec!["-s", "S1", "-H", "R"],
                SampleMode::SingleSample {
                    idx: 0,
                    spec: HaplotypeSpec::parse("R").unwrap(),
                },
            ),
            (
                vec!["-s", "S1", "-H", "A"],
                SampleMode::SingleSample {
                    idx: 0,
                    spec: HaplotypeSpec::parse("A").unwrap(),
                },
            ),
            (
                vec!["-s", "S1", "-H", "I"],
                SampleMode::SingleSample {
                    idx: 0,
                    spec: HaplotypeSpec::parse("I").unwrap(),
                },
            ),
            (
                vec!["-s", "S1", "-H", "1"],
                SampleMode::SingleSample {
                    idx: 0,
                    spec: HaplotypeSpec::parse("1").unwrap(),
                },
            ),
            (
                vec!["-s", "S1", "-H", "2"],
                SampleMode::SingleSample {
                    idx: 0,
                    spec: HaplotypeSpec::parse("2").unwrap(),
                },
            ),
            (
                vec!["-s", "S1", "-H", "1pIu"],
                SampleMode::SingleSample {
                    idx: 0,
                    spec: HaplotypeSpec::parse("1pIu").unwrap(),
                },
            ),
            (
                vec!["-s", "S1", "-H", "2pIu"],
                SampleMode::SingleSample {
                    idx: 0,
                    spec: HaplotypeSpec::parse("2pIu").unwrap(),
                },
            ),
            (
                vec!["-s", "S1", "-H", "LR"],
                SampleMode::SingleSample {
                    idx: 0,
                    spec: HaplotypeSpec::parse("LR").unwrap(),
                },
            ),
            (
                vec!["-s", "S1", "-H", "LA"],
                SampleMode::SingleSample {
                    idx: 0,
                    spec: HaplotypeSpec::parse("LA").unwrap(),
                },
            ),
            (
                vec!["-s", "S1", "-H", "SR"],
                SampleMode::SingleSample {
                    idx: 0,
                    spec: HaplotypeSpec::parse("SR").unwrap(),
                },
            ),
            (
                vec!["-s", "S1", "-H", "SA"],
                SampleMode::SingleSample {
                    idx: 0,
                    spec: HaplotypeSpec::parse("SA").unwrap(),
                },
            ),
            // -s S2 (no -H): IUPAC across the one sample
            (
                vec!["-s", "S2"],
                SampleMode::IupacAllSamples { samples: vec![1] },
            ),
        ];

        for (cli, mode) in &modes {
            let mut args = vec!["consensus", "-f"];
            args.push(ref_fa.to_str().unwrap());
            args.extend(cli.iter());
            args.push(gz.to_str().unwrap());
            let out = Command::new("bcftools")
                .args(&args)
                .output()
                .expect("bcftools");
            assert!(
                out.status.success(),
                "bcftools {:?} failed: {}",
                cli,
                String::from_utf8_lossy(&out.stderr)
            );
            let theirs: Vec<u8> = String::from_utf8_lossy(&out.stdout)
                .lines()
                .filter(|l| !l.starts_with('>') && !l.is_empty())
                .flat_map(|l| l.as_bytes().iter().copied())
                .collect();

            let opts = ApplyOptions {
                sample_mode: mode.clone(),
                ..Default::default()
            };
            let mut state = ApplyState::new(0, seq.as_bytes().to_vec());
            for r in &recs {
                apply_variant(r, &mut state, &opts, None);
            }
            assert_eq!(
                state.buf,
                theirs,
                "cli {:?} mismatch:\n  ours:   {}\n  theirs: {}",
                cli,
                String::from_utf8_lossy(&state.buf),
                String::from_utf8_lossy(&theirs)
            );
        }

        // ialt==0 overlap REF trim (consensus.c:811-825): a deletion followed by
        // a hom-ref variant whose REF spans the frozen position. With -s S1 (no
        // -H), hom-ref GT → ialt=0, which must trim forward and apply REF.
        {
            let body = format!(
                "##fileformat=VCFv4.3\n\
                 ##contig=<ID=chr1,length=100>\n\
                 ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"G\">\n\
                 #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
                 chr1\t4\t.\t{r4}{r5}{r6}\t{r4}\t.\t.\t.\tGT\t0/1\n\
                 chr1\t5\t.\t{r5}{r6}{r7}\t{r5}\t.\t.\t.\tGT\t0/0\n",
                r4 = r(4),
                r5 = r(5),
                r6 = r(6),
                r7 = r(7),
            );
            let vcf = dir.join("ialt0.vcf");
            std::fs::write(&vcf, body).unwrap();
            let st = Command::new("bgzip").arg(&vcf).status().expect("bgzip");
            assert!(st.success());
            let gz2 = vcf.with_extension("vcf.gz");
            let st = Command::new("tabix")
                .arg("-p")
                .arg("vcf")
                .arg(&gz2)
                .status()
                .expect("tabix");
            assert!(st.success());

            let out = Command::new("bcftools")
                .args([
                    "consensus",
                    "-f",
                    ref_fa.to_str().unwrap(),
                    "-s",
                    "S1",
                    gz2.to_str().unwrap(),
                ])
                .output()
                .expect("bcftools");
            assert!(
                out.status.success(),
                "ialt0 bcftools failed: {}",
                String::from_utf8_lossy(&out.stderr)
            );
            let theirs: Vec<u8> = String::from_utf8_lossy(&out.stdout)
                .lines()
                .filter(|l| !l.starts_with('>') && !l.is_empty())
                .flat_map(|l| l.as_bytes().iter().copied())
                .collect();

            let store2 = VcfStore::load(&gz2).unwrap();
            let recs2 = store2.query("chr1", 0, end, 1);
            let s1 = store2.sample_index("S1").unwrap();
            let opts = ApplyOptions {
                sample_mode: SampleMode::IupacAllSamples { samples: vec![s1] },
                ..Default::default()
            };
            let mut state = ApplyState::new(0, seq.as_bytes().to_vec());
            for r in &recs2 {
                apply_variant(r, &mut state, &opts, None);
            }
            assert_eq!(
                state.buf,
                theirs,
                "ialt0 trim mismatch:\n  ours:   {}\n  theirs: {}",
                String::from_utf8_lossy(&state.buf),
                String::from_utf8_lossy(&theirs)
            );
        }
    }

    /// M4 acceptance: parity for `-a`/`-M`/`--mark-*`/`-m`/`-c`. Ignored by
    /// default; run with `cargo test -- --ignored bcftools_m4_parity`.
    #[test]
    #[ignore]
    fn bcftools_m4_parity() {
        use crate::chain::Chain;
        use crate::mask::{Mask, MaskWith};
        use std::process::Command;
        let dir = std::env::temp_dir().join("consensus_rs_m4_parity");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let seqb = seq.as_bytes();
        std::fs::write(dir.join("ref.fa"), format!(">chr1\n{}\n", seq)).unwrap();
        let ref_fa = dir.join("ref.fa");
        let r = |p1: usize| -> char { seqb[p1 - 1] as char };

        // One sample, several variant types + a missing GT, with a gap for -a.
        let body = format!(
            "##fileformat=VCFv4.3\n\
             ##contig=<ID=chr1,length=100>\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"G\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
             chr1\t2\t.\t{r2}\tG\t.\t.\t.\tGT\t0/1\n\
             chr1\t10\t.\t{r10}\t{r10}TTT\t.\t.\t.\tGT\t0/1\n\
             chr1\t15\t.\t{r15}{r16}\t{r15}\t.\t.\t.\tGT\t0/1\n\
             chr1\t40\t.\t{r40}\tA\t.\t.\t.\tGT\t./.\n",
            r2 = r(2),
            r10 = r(10),
            r15 = r(15),
            r16 = r(16),
            r40 = r(40),
        );
        let vcf = dir.join("m4.vcf");
        std::fs::write(&vcf, body).unwrap();
        let st = Command::new("bgzip").arg(&vcf).status().expect("bgzip");
        assert!(st.success());
        let gz = vcf.with_extension("vcf.gz");
        let st = Command::new("tabix")
            .arg("-p")
            .arg("vcf")
            .arg(&gz)
            .status()
            .expect("tabix");
        assert!(st.success());

        let store = VcfStore::load(&gz).unwrap();
        let end = seq.len() as i64 - 1;
        let recs = store.query("chr1", 0, end, 1);
        let s1 = store.sample_index("S1").unwrap();

        // mask bed: chr1 30-32 (0-based [30,33)) — sits after the ins(10)/del(15)
        // so mod_off!=0 when the mask is applied, exercising the coord fix.
        std::fs::write(dir.join("mask.bed"), "chr1\t30\t33\n").unwrap();
        let mask_bed = dir.join("mask.bed");

        // helper to run bcftools (on the main gz) and strip header
        let bcftools_seq = |cli: &[&str]| -> Vec<u8> {
            let mut args = vec!["consensus", "-f"];
            args.push(ref_fa.to_str().unwrap());
            args.extend(cli.iter());
            args.push(gz.to_str().unwrap());
            let out = Command::new("bcftools")
                .args(&args)
                .output()
                .expect("bcftools");
            assert!(
                out.status.success(),
                "bcftools {:?} failed: {}",
                cli,
                String::from_utf8_lossy(&out.stderr)
            );
            String::from_utf8_lossy(&out.stdout)
                .lines()
                .filter(|l| !l.starts_with('>') && !l.is_empty())
                .flat_map(|l| l.as_bytes().iter().copied())
                .collect()
        };
        // (cli, opts_builder). Each builds ApplyOptions + optional mask/chain.
        // 1. -a N (absent)
        {
            let theirs = bcftools_seq(&["-a", "N"]);
            let opts = ApplyOptions {
                absent_allele: Some(b'N'),
                sample_mode: SampleMode::IupacAllSamples { samples: vec![s1] },
                ..Default::default()
            };
            let mut state = ApplyState::new(0, seq.as_bytes().to_vec()).with_chr("chr1");
            for r in &recs {
                apply_variant(r, &mut state, &opts, None);
            }
            if let Some(a) = opts.absent_allele {
                apply_absent(&mut state, i64::MAX, a);
            }
            assert_eq!(state.buf, theirs, "-a N");
        }
        // 2. -M ? (missing char) with single sample -H 1 (so missing becomes ?)
        {
            let theirs = bcftools_seq(&["-s", "S1", "-H", "1", "-M", "?"]);
            let opts = ApplyOptions {
                missing_allele: Some(b'?'),
                sample_mode: SampleMode::SingleSample {
                    idx: s1,
                    spec: crate::haplotype::HaplotypeSpec::parse("1").unwrap(),
                },
                ..Default::default()
            };
            let mut state = ApplyState::new(0, seq.as_bytes().to_vec()).with_chr("chr1");
            for r in &recs {
                apply_variant(r, &mut state, &opts, None);
            }
            assert_eq!(state.buf, theirs, "-M ?");
        }
        // 3. --mark-ins uc
        {
            let theirs = bcftools_seq(&["--mark-ins", "uc"]);
            let opts = ApplyOptions {
                mark_ins: Some(TO_UPPER as u8),
                sample_mode: SampleMode::IupacAllSamples { samples: vec![s1] },
                ..Default::default()
            };
            let mut state = ApplyState::new(0, seq.as_bytes().to_vec()).with_chr("chr1");
            for r in &recs {
                apply_variant(r, &mut state, &opts, None);
            }
            assert_eq!(state.buf, theirs, "--mark-ins uc");
        }
        // 4. --mark-snv lc
        {
            let theirs = bcftools_seq(&["--mark-snv", "lc"]);
            let opts = ApplyOptions {
                mark_snv: Some(TO_LOWER as u8),
                sample_mode: SampleMode::IupacAllSamples { samples: vec![s1] },
                ..Default::default()
            };
            let mut state = ApplyState::new(0, seq.as_bytes().to_vec()).with_chr("chr1");
            for r in &recs {
                apply_variant(r, &mut state, &opts, None);
            }
            assert_eq!(state.buf, theirs, "--mark-snv lc");
        }
        // 5. --mark-del #
        {
            let theirs = bcftools_seq(&["--mark-del", "#"]);
            let opts = ApplyOptions {
                mark_del: Some(b'#'),
                sample_mode: SampleMode::IupacAllSamples { samples: vec![s1] },
                ..Default::default()
            };
            let mut state = ApplyState::new(0, seq.as_bytes().to_vec()).with_chr("chr1");
            for r in &recs {
                apply_variant(r, &mut state, &opts, None);
            }
            assert_eq!(state.buf, theirs, "--mark-del #");
        }
        // 6. -m mask.bed (char N, skips variants). Uses the indel-bearing main
        // VCF: mask bed [30,33) sits after the ins(10)/del(15), so mod_off!=0 —
        // this is the case the old post-apply mask mislocated.
        {
            let theirs = bcftools_seq(&["-m", mask_bed.to_str().unwrap()]);
            let opts = ApplyOptions {
                mask: Some(Rc::new(
                    Mask::load(&mask_bed, MaskWith::Char(b'N')).unwrap(),
                )),
                sample_mode: SampleMode::IupacAllSamples { samples: vec![s1] },
                ..Default::default()
            };
            let state = apply_region("chr1", seq.as_bytes().to_vec(), 0, &recs, &opts, None);
            assert_eq!(state.buf, theirs, "-m mask N (indel before mask)");
        }
        // 7. -m mask.bed --mask-with lc (no variant skip) — indel VCF
        {
            let theirs = bcftools_seq(&["-m", mask_bed.to_str().unwrap(), "--mask-with", "lc"]);
            let opts = ApplyOptions {
                mask: Some(Rc::new(Mask::load(&mask_bed, MaskWith::Lc).unwrap())),
                sample_mode: SampleMode::IupacAllSamples { samples: vec![s1] },
                ..Default::default()
            };
            let state = apply_region("chr1", seq.as_bytes().to_vec(), 0, &recs, &opts, None);
            assert_eq!(state.buf, theirs, "-m mask lc (indel before mask)");
        }
        // 8. -c chain (compare chain text)
        {
            let chain_path = dir.join("ours.chain");
            let theirs_chain = {
                let p = dir.join("theirs.chain");
                let out = Command::new("bcftools")
                    .args([
                        "consensus",
                        "-f",
                        ref_fa.to_str().unwrap(),
                        "-c",
                        p.to_str().unwrap(),
                        gz.to_str().unwrap(),
                    ])
                    .output()
                    .expect("bcftools");
                assert!(out.status.success());
                std::fs::read_to_string(&p).unwrap()
            };
            let opts = ApplyOptions {
                sample_mode: SampleMode::IupacAllSamples { samples: vec![s1] },
                ..Default::default()
            };
            let mut chain = Chain::new("chr1".into(), 0, seq.len() as i64);
            let _state = apply_region(
                "chr1",
                seq.as_bytes().to_vec(),
                0,
                &recs,
                &opts,
                Some(&mut chain),
            );
            let ours = chain.render();
            std::fs::write(&chain_path, &ours).unwrap();
            assert_eq!(ours, theirs_chain, "-c chain");
        }
    }
}
