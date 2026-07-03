//! engine — ConsensusEngine: holds preprocessed material (ref + vcfs), runs
//! multi-threaded consensus production.
//!
//! (docs/design.md §5.4 / §5.6) Each (region, sample, haplotype) task is
//! independent → natural parallelism. `consensus_many` and `consensus_iter`
//! group identical regions so ref fetch, VCF query, and region planning are
//! amortized across sample/haplotype tasks.
//!
//! This module is PyO3-free; `py.rs` wraps it under the `python` feature.

use crate::apply::{apply_region_planned_set, ApplyOptions, TO_LOWER, TO_UPPER};
use crate::chain::Chain;
use crate::compiled::{
    allele_case_flags, RecordFlags, ALLELE_HAS_ASCII_LOWER, ALLELE_HAS_ASCII_UPPER,
};
use crate::haplotype::{HaplotypeSpec, SampleMode};
use crate::planner::{plan_region_set, PlanOptions, RegionPlan};
use crate::ref_index::RefIndex;
use crate::stats::FastPathLane;
use crate::vcf_store::{LoadStrategy, RecordSet, VcfStore};
use crossbeam_channel::{bounded, Receiver, Sender};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet, VecDeque};
use std::path::PathBuf;
use std::sync::Arc;

/// One production task: a region + how to pick alleles for it.
#[derive(Clone)]
pub struct ConsensusTask {
    pub chr: String,
    /// 1-based inclusive start.
    pub start: i64,
    /// 1-based inclusive end.
    pub end: i64,
    pub vcf_key: String,
    pub gene_id: String,
    pub sample: Option<String>,
    pub haplotype: Option<String>,
}

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
struct TaskGroupKey {
    chr: String,
    start: i64,
    end: i64,
    vcf_key: String,
}

struct TaskGroup {
    key: TaskGroupKey,
    indices: Vec<usize>,
}

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
struct TaskExecKey {
    sample: Option<String>,
    haplotype: Option<String>,
}

#[derive(Clone)]
struct CachedOutput {
    seq: Vec<u8>,
    chain: Option<String>,
    error: Option<String>,
}

#[derive(Clone, Copy)]
struct BatchTask {
    sample_idx: usize,
    hap: usize,
}

struct RefOnlyBatchPatch<'a> {
    idx: usize,
    rlen: usize,
    ref_allele: &'a [u8],
    ref_case_flags: u8,
    ref_out: u8,
    to_upper: bool,
}

struct Snp1BatchPatch<'a> {
    idx: usize,
    ref_out: u8,
    alt_out: u8,
    gt_bits: &'a crate::vcf_store::BiallelicPhasedGtBits,
}

struct MnpBatchPatch<'a> {
    idx: usize,
    rlen: usize,
    ref_allele: &'a [u8],
    ref_case_flags: u8,
    alt: &'a [u8],
    alt_case_flags: u8,
    gt_bits: &'a crate::vcf_store::BiallelicPhasedGtBits,
    to_upper: bool,
}

enum SameLenBatchPatch<'a> {
    RefOnly(RefOnlyBatchPatch<'a>),
    Snp1(Snp1BatchPatch<'a>),
    Mnp(MnpBatchPatch<'a>),
}

/// One produced result.
pub struct ConsensusResult {
    pub gene_id: String,
    pub sample: Option<String>,
    pub haplotype: Option<String>,
    pub seq: Vec<u8>,
    pub chain: Option<String>,
    pub error: Option<String>,
}

/// Engine config (cli-style options applied to every task).
#[derive(Clone, Default)]
pub struct EngineOptions {
    pub iupac_codes: bool,
    pub missing: Option<u8>,
    pub absent: Option<u8>,
    pub mark_del: Option<u8>,
    pub mark_ins: Option<u8>,
    pub mark_snv: Option<u8>,
    pub mask: Option<PathBuf>,
    pub mask_with: crate::mask::MaskWith,
    pub chain: bool,
    pub regions_overlap: u8,
}

/// The engine: preprocessed ref + a map of VcfStores, Send+Sync via Arc.
/// Cheap to clone (all fields are Arc / Clone).
#[derive(Clone)]
pub struct ConsensusEngine {
    ref_index: Arc<RefIndex>,
    vcfs: Arc<HashMap<String, Arc<VcfStore>>>,
    mask: Option<Arc<crate::mask::Mask>>,
    opts: EngineOptions,
}

impl ConsensusEngine {
    pub fn new(ref_index: RefIndex, vcfs: HashMap<String, VcfStore>, opts: EngineOptions) -> Self {
        Self::try_new(ref_index, vcfs, opts).expect("failed to load mask")
    }

    pub fn try_new(
        ref_index: RefIndex,
        vcfs: HashMap<String, VcfStore>,
        opts: EngineOptions,
    ) -> Result<Self, String> {
        let vcfs = vcfs.into_iter().map(|(k, v)| (k, Arc::new(v))).collect();
        let mask = match &opts.mask {
            Some(path) => Some(Arc::new(crate::mask::Mask::load(path, opts.mask_with)?)),
            None => None,
        };
        Ok(ConsensusEngine {
            ref_index: Arc::new(ref_index),
            vcfs: Arc::new(vcfs),
            mask,
            opts,
        })
    }

    /// Eagerly load ref + all VCFs.
    pub fn load(
        ref_path: impl Into<PathBuf>,
        vcf_paths: HashMap<String, PathBuf>,
        opts: EngineOptions,
    ) -> Result<Self, String> {
        let ref_index = RefIndex::load(ref_path)?;
        let mut vcfs = HashMap::new();
        for (k, p) in vcf_paths {
            vcfs.insert(k, VcfStore::load_with_strategy(p, LoadStrategy::Eager)?);
        }
        ConsensusEngine::try_new(ref_index, vcfs, opts)
    }

    /// Run all tasks in parallel, returning results in input order.
    ///
    /// `threads` is the size of the rayon pool created for this call. The
    /// engine holds no compute resources of its own — only the preprocessed
    /// ref + VCFs — so the pool is built and torn down per call, and the
    /// thread count is the caller's to decide each time.
    pub fn consensus_many(
        &self,
        tasks: Vec<ConsensusTask>,
        threads: usize,
    ) -> Vec<ConsensusResult> {
        let n = tasks.len();
        if n == 0 {
            return Vec::new();
        }
        let nthr = threads.max(1);
        let pool = thread_pool(nthr);

        let groups = group_tasks(&tasks);
        let tasks = Arc::new(tasks);
        let engine = self.clone_shallow();
        let indexed: Vec<Vec<(usize, ConsensusResult)>> = pool.install(move || {
            groups
                .par_iter()
                .map(|group| engine.run_group(&tasks, group))
                .collect()
        });

        let mut out: Vec<Option<ConsensusResult>> = (0..n).map(|_| None).collect();
        for group_results in indexed {
            for (idx, result) in group_results {
                out[idx] = Some(result);
            }
        }
        out.into_iter().flatten().collect()
    }

    /// Run a single task (used both standalone and inside consensus_many).
    pub fn run_one(&self, task: &ConsensusTask) -> ConsensusResult {
        // Build an error result by cloning task labels (Fn, not FnOnce).
        let mk_err = |task: &ConsensusTask, err: String| -> ConsensusResult {
            ConsensusResult {
                gene_id: task.gene_id.clone(),
                sample: task.sample.clone(),
                haplotype: task.haplotype.clone(),
                seq: Vec::new(),
                chain: None,
                error: Some(err),
            }
        };

        let vcf = match self.vcfs.get(&task.vcf_key) {
            Some(v) => v,
            None => return mk_err(task, format!("unknown vcf_key: {}", task.vcf_key)),
        };

        // Fetch ref [start, end] (1-based inclusive) -> plus-strand bytes.
        let ref_seq = match self.ref_index.fetch_1based(&task.chr, task.start, task.end) {
            Ok(s) => s,
            Err(e) => return mk_err(task, format!("ref fetch failed: {}", e)),
        };
        let ori_pos = task.start - 1; // 0-based
        let end0 = task.end - 1;

        // Region query + fastpath plan (0-based inclusive). The plan is a
        // conservative dispatch hint; each fastpath still validates before it
        // mutates output.
        let records = vcf.query_set(&task.chr, ori_pos, end0, self.opts.regions_overlap);
        let plan_opts = self.plan_options_for_records(&task.chr, &records);
        let plan = plan_region_set(&records, plan_opts);

        // Build sample mode + options for this task.
        let sample_mode =
            build_sample_mode(vcf, &task.sample, &task.haplotype, self.opts.iupac_codes);
        let opts = ApplyOptions {
            absent_allele: self.opts.absent,
            missing_allele: self.opts.missing,
            mark_del: self.opts.mark_del,
            mark_ins: self.opts.mark_ins,
            mark_snv: self.opts.mark_snv,
            sample_mode,
            mask: self.mask.clone(),
        };

        if self.opts.chain {
            let mut chain = Chain::new(task.chr.clone(), ori_pos, ref_seq.len() as i64);
            let state = apply_region_planned_set(
                &task.chr,
                ref_seq,
                ori_pos,
                &records,
                &opts,
                Some(&mut chain),
                Some(&plan),
            );
            ConsensusResult {
                gene_id: task.gene_id.clone(),
                sample: task.sample.clone(),
                haplotype: task.haplotype.clone(),
                seq: state.buf,
                chain: Some(chain.render()),
                error: None,
            }
        } else {
            let state = apply_region_planned_set(
                &task.chr,
                ref_seq,
                ori_pos,
                &records,
                &opts,
                None,
                Some(&plan),
            );
            ConsensusResult {
                gene_id: task.gene_id.clone(),
                sample: task.sample.clone(),
                haplotype: task.haplotype.clone(),
                seq: state.buf,
                chain: None,
                error: None,
            }
        }
    }

    fn run_group(
        &self,
        tasks: &[ConsensusTask],
        group: &TaskGroup,
    ) -> Vec<(usize, ConsensusResult)> {
        let first_idx = group.indices[0];
        let first_task = &tasks[first_idx];
        let group_error = |err: String| -> Vec<(usize, ConsensusResult)> {
            group
                .indices
                .iter()
                .map(|&i| (i, error_result(&tasks[i], err.clone())))
                .collect()
        };

        let vcf = match self.vcfs.get(&group.key.vcf_key) {
            Some(v) => v,
            None => return group_error(format!("unknown vcf_key: {}", group.key.vcf_key)),
        };

        let ref_seq =
            match self
                .ref_index
                .fetch_1based(&group.key.chr, group.key.start, group.key.end)
            {
                Ok(s) => s,
                Err(e) => return group_error(format!("ref fetch failed: {}", e)),
            };
        let ori_pos = group.key.start - 1;
        let end0 = group.key.end - 1;
        let records = vcf.query_set(&group.key.chr, ori_pos, end0, self.opts.regions_overlap);
        let plan_opts = self.plan_options_for_records(&group.key.chr, &records);
        let plan = plan_region_set(&records, plan_opts);
        if let Some(batch) = self
            .try_run_biallelic_phased_batch(tasks, group, vcf, &ref_seq, ori_pos, &records, &plan)
        {
            return batch;
        }
        let shared_mask = self.mask.clone();
        if records.is_empty() {
            return self.run_empty_group(tasks, group, &ref_seq, ori_pos, shared_mask.as_deref());
        }

        let cache_enabled = has_duplicate_exec_keys(tasks, &group.indices);
        let mut output_cache: HashMap<TaskExecKey, CachedOutput> = HashMap::new();
        let mut ref_seq = Some(ref_seq);
        let n = group.indices.len();
        let mut out = Vec::with_capacity(n);
        for (j, &idx) in group.indices.iter().enumerate() {
            let task = &tasks[idx];
            debug_assert_eq!(task.chr, first_task.chr);
            debug_assert_eq!(task.start, first_task.start);
            debug_assert_eq!(task.end, first_task.end);
            debug_assert_eq!(task.vcf_key, first_task.vcf_key);

            let cache_key = if cache_enabled {
                let key = task_exec_key(task);
                if let Some(cached) = output_cache.get(&key) {
                    out.push((idx, result_from_cached(task, cached)));
                    continue;
                }
                Some(key)
            } else {
                None
            };

            let ref_for_task = if cache_enabled {
                ref_seq.as_ref().expect("shared ref available").clone()
            } else if j + 1 == n {
                ref_seq.take().expect("last task consumes shared ref")
            } else {
                ref_seq.as_ref().expect("shared ref available").clone()
            };

            let result = self.run_group_task(
                task,
                vcf,
                ref_for_task,
                ori_pos,
                &records,
                &plan,
                shared_mask.clone(),
            );
            if let Some(key) = cache_key {
                output_cache.insert(key, CachedOutput::from(&result));
            }
            out.push((idx, result));
        }
        out
    }

    fn run_empty_group(
        &self,
        tasks: &[ConsensusTask],
        group: &TaskGroup,
        ref_seq: &[u8],
        ori_pos: i64,
        shared_mask: Option<&crate::mask::Mask>,
    ) -> Vec<(usize, ConsensusResult)> {
        let mut seq = match self.opts.absent {
            Some(absent) => vec![absent; ref_seq.len()],
            None => ref_seq.to_vec(),
        };
        if self.opts.absent.is_none() {
            if let Some(mask) = shared_mask {
                mask.apply_to_buf(&group.key.chr, &mut seq, ori_pos);
            }
        }
        let chain = if self.opts.chain {
            let mut chain = Chain::new(group.key.chr.clone(), ori_pos, ref_seq.len() as i64);
            Some(chain.render())
        } else {
            None
        };

        let mut seq = Some(seq);
        let mut out = Vec::with_capacity(group.indices.len());
        let n = group.indices.len();
        for (j, &idx) in group.indices.iter().enumerate() {
            let task = &tasks[idx];
            let seq_out = if j + 1 == n {
                seq.take().expect("last empty-group result consumes seq")
            } else {
                seq.as_ref()
                    .expect("empty-group template seq available")
                    .clone()
            };
            out.push((
                idx,
                ConsensusResult {
                    gene_id: task.gene_id.clone(),
                    sample: task.sample.clone(),
                    haplotype: task.haplotype.clone(),
                    seq: seq_out,
                    chain: chain.clone(),
                    error: None,
                },
            ));
        }
        out
    }

    #[allow(clippy::too_many_arguments)]
    fn run_group_task(
        &self,
        task: &ConsensusTask,
        vcf: &VcfStore,
        ref_for_task: Vec<u8>,
        ori_pos: i64,
        records: &RecordSet<'_>,
        plan: &RegionPlan,
        shared_mask: Option<Arc<crate::mask::Mask>>,
    ) -> ConsensusResult {
        let sample_mode =
            build_sample_mode(vcf, &task.sample, &task.haplotype, self.opts.iupac_codes);
        let opts = ApplyOptions {
            absent_allele: self.opts.absent,
            missing_allele: self.opts.missing,
            mark_del: self.opts.mark_del,
            mark_ins: self.opts.mark_ins,
            mark_snv: self.opts.mark_snv,
            sample_mode,
            mask: shared_mask,
        };

        if self.opts.chain {
            let mut chain = Chain::new(task.chr.clone(), ori_pos, ref_for_task.len() as i64);
            let state = apply_region_planned_set(
                &task.chr,
                ref_for_task,
                ori_pos,
                records,
                &opts,
                Some(&mut chain),
                Some(plan),
            );
            ConsensusResult {
                gene_id: task.gene_id.clone(),
                sample: task.sample.clone(),
                haplotype: task.haplotype.clone(),
                seq: state.buf,
                chain: Some(chain.render()),
                error: None,
            }
        } else {
            let state = apply_region_planned_set(
                &task.chr,
                ref_for_task,
                ori_pos,
                records,
                &opts,
                None,
                Some(plan),
            );
            ConsensusResult {
                gene_id: task.gene_id.clone(),
                sample: task.sample.clone(),
                haplotype: task.haplotype.clone(),
                seq: state.buf,
                chain: None,
                error: None,
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn try_run_biallelic_phased_batch(
        &self,
        tasks: &[ConsensusTask],
        group: &TaskGroup,
        vcf: &VcfStore,
        ref_seq: &[u8],
        ori_pos: i64,
        records: &RecordSet<'_>,
        plan: &RegionPlan,
    ) -> Option<Vec<(usize, ConsensusResult)>> {
        if records.is_empty() || plan.lane != FastPathLane::SameLenOnly || self.opts.chain {
            return None;
        }
        let base_ref_storage = self.mask.as_ref().map(|mask| {
            let mut masked = ref_seq.to_vec();
            mask.apply_to_buf(&group.key.chr, &mut masked, ori_pos);
            masked
        });
        let base_ref = base_ref_storage.as_deref().unwrap_or(ref_seq);

        let n_samples = vcf.n_sample() as usize;
        let n_words = n_samples.div_ceil(64);
        let mut task_by_hap_sample = [vec![usize::MAX; n_samples], vec![usize::MAX; n_samples]];
        let mut active_words_by_hap = [vec![0u64; n_words], vec![0u64; n_words]];
        let mut active_word_indices_by_hap = [Vec::new(), Vec::new()];
        let mut batch_tasks = Vec::with_capacity(group.indices.len());
        let mut output_map = Vec::with_capacity(group.indices.len());
        for &input_idx in &group.indices {
            let task = &tasks[input_idx];
            let sample = task.sample.as_ref()?;
            let sample_idx = vcf.sample_index(sample)?;
            if sample_idx < 0 {
                return None;
            }
            let sample_idx = sample_idx as usize;
            if sample_idx >= n_samples {
                return None;
            }
            let hap = task
                .haplotype
                .as_ref()
                .and_then(|h| parse_haplotype_index_no_alloc(h))?;
            if hap == 0 || hap > 2 {
                return None;
            }
            let hap_idx = hap as usize - 1;
            let word_idx = sample_idx / 64;
            if active_words_by_hap[hap_idx][word_idx] == 0 {
                active_word_indices_by_hap[hap_idx].push(word_idx);
            }
            active_words_by_hap[hap_idx][word_idx] |= 1u64 << (sample_idx & 63);
            let slot = &mut task_by_hap_sample[hap_idx][sample_idx];
            let batch_idx = if *slot == usize::MAX {
                let idx = batch_tasks.len();
                *slot = idx;
                batch_tasks.push(BatchTask {
                    sample_idx,
                    hap: hap as usize,
                });
                idx
            } else {
                *slot
            };
            output_map.push((input_idx, batch_idx));
        }
        if batch_tasks.is_empty() {
            return None;
        }

        let mut patches = Vec::with_capacity(records.len());
        let mut frz_pos = -1i64;
        for rec in records.iter() {
            if rec.alleles.len() == 1 {
                if self.opts.absent.is_none() {
                    continue;
                }
                if rec.pos <= frz_pos || rec.rlen <= 0 {
                    return None;
                }
                let rlen = rec.rlen as usize;
                if rec.alleles[0].len() != rlen {
                    return None;
                }
                let idx = rec.pos - ori_pos;
                if idx < 0 {
                    return None;
                }
                let idx = idx as usize;
                if idx + rlen > base_ref.len() {
                    return None;
                }
                if !base_ref[idx..idx + rlen].eq_ignore_ascii_case(&rec.alleles[0]) {
                    return None;
                }
                let ref_case_flags = rec
                    .compiled
                    .allele_op(0)
                    .map(|op| op.case_flags)
                    .unwrap_or_else(|| allele_case_flags(&rec.alleles[0]));
                let to_upper = base_ref[idx].is_ascii_uppercase();
                patches.push(SameLenBatchPatch::RefOnly(RefOnlyBatchPatch {
                    idx,
                    rlen,
                    ref_allele: &rec.alleles[0],
                    ref_case_flags,
                    ref_out: snp1_alt_with_case_and_mark(
                        rec.alleles[0][0],
                        rec.alleles[0][0],
                        to_upper,
                        ref_case_flags,
                        None,
                    ),
                    to_upper,
                }));
                frz_pos = rec.ref_end();
                continue;
            }
            if rec.pos <= frz_pos
                || rec.rlen <= 0
                || !rec.compiled.flags.contains(RecordFlags::BIALLELIC)
                || !rec.compiled.flags.contains(RecordFlags::ALL_ALT_SAME_LEN)
            {
                return None;
            }
            let gt_bits = rec.gt_bits.as_ref()?;
            let rlen = rec.rlen as usize;
            if rec.alleles.len() != 2
                || rec.alleles[0].len() != rlen
                || rec.alleles[1].len() != rlen
            {
                return None;
            }
            let idx = rec.pos - ori_pos;
            if idx < 0 {
                return None;
            }
            let idx = idx as usize;
            if idx + rlen > base_ref.len() {
                return None;
            }
            if !base_ref[idx..idx + rlen].eq_ignore_ascii_case(&rec.alleles[0]) {
                return None;
            }
            let ref_case_flags = rec
                .compiled
                .allele_op(0)
                .map(|op| op.case_flags)
                .unwrap_or_else(|| allele_case_flags(&rec.alleles[0]));
            let alt_case_flags = rec
                .compiled
                .allele_op(1)
                .map(|op| op.case_flags)
                .unwrap_or_else(|| allele_case_flags(&rec.alleles[1]));
            let to_upper = base_ref[idx].is_ascii_uppercase();
            if rlen == 1 {
                let ref_base = rec.alleles[0][0];
                let alt_base = rec.alleles[1][0];
                patches.push(SameLenBatchPatch::Snp1(Snp1BatchPatch {
                    idx,
                    ref_out: snp1_alt_with_case_and_mark(
                        ref_base,
                        ref_base,
                        to_upper,
                        ref_case_flags,
                        None,
                    ),
                    alt_out: snp1_alt_with_case_and_mark(
                        ref_base,
                        alt_base,
                        to_upper,
                        alt_case_flags,
                        self.opts.mark_snv,
                    ),
                    gt_bits,
                }));
            } else {
                patches.push(SameLenBatchPatch::Mnp(MnpBatchPatch {
                    idx,
                    rlen,
                    ref_allele: &rec.alleles[0],
                    ref_case_flags,
                    alt: &rec.alleles[1],
                    alt_case_flags,
                    gt_bits,
                    to_upper,
                }));
            }
            frz_pos = rec.ref_end();
        }
        if patches.is_empty() {
            return None;
        }

        if let Some(out) = self.try_run_biallelic_phased_alt_only_batch(
            tasks,
            base_ref,
            &output_map,
            &task_by_hap_sample,
            &active_words_by_hap,
            &active_word_indices_by_hap,
            &patches,
        ) {
            return Some(out);
        }
        if let Some(out) = self.try_run_biallelic_phased_missing_batch(
            tasks,
            base_ref,
            &output_map,
            &task_by_hap_sample,
            &active_words_by_hap,
            &active_word_indices_by_hap,
            &patches,
        ) {
            return Some(out);
        }
        if let Some(out) = self.try_run_biallelic_phased_absent_batch(
            tasks,
            base_ref,
            &output_map,
            &task_by_hap_sample,
            &active_words_by_hap,
            &active_word_indices_by_hap,
            &patches,
        ) {
            return Some(out);
        }

        let mut buffers: Vec<Vec<u8>> = match self.opts.absent {
            Some(absent) => batch_tasks
                .iter()
                .map(|_| vec![absent; base_ref.len()])
                .collect(),
            None => batch_tasks.iter().map(|_| base_ref.to_vec()).collect(),
        };
        for patch in &patches {
            match patch {
                SameLenBatchPatch::RefOnly(patch) => {
                    if self.opts.absent.is_none() {
                        continue;
                    }
                    for buf in &mut buffers {
                        if patch.rlen == 1 {
                            buf[patch.idx] = patch.ref_out;
                        } else {
                            let dst = &mut buf[patch.idx..patch.idx + patch.rlen];
                            copy_alt_with_case_flags(
                                dst,
                                patch.ref_allele,
                                patch.to_upper,
                                patch.ref_case_flags,
                            );
                        }
                    }
                }
                SameLenBatchPatch::Snp1(patch) => {
                    for (task_i, task) in batch_tasks.iter().enumerate() {
                        let buf = &mut buffers[task_i];
                        match patch.gt_bits.allele_for_hap(task.sample_idx, task.hap)? {
                            Some(0) if self.opts.absent.is_some() => {
                                buf[patch.idx] = patch.ref_out;
                            }
                            Some(0) => {}
                            Some(1) => {
                                buf[patch.idx] = patch.alt_out;
                            }
                            None => {
                                if let Some(missing) = self.opts.missing {
                                    buf[patch.idx] = missing;
                                }
                            }
                            _ => return None,
                        }
                    }
                }
                SameLenBatchPatch::Mnp(patch) => {
                    for (task_i, task) in batch_tasks.iter().enumerate() {
                        let allele =
                            match patch.gt_bits.allele_for_hap(task.sample_idx, task.hap)? {
                                Some(0) if self.opts.absent.is_some() => {
                                    Some(Ok((patch.ref_allele, patch.ref_case_flags)))
                                }
                                Some(0) => None,
                                Some(1) => Some(Ok((patch.alt, patch.alt_case_flags))),
                                None => self.opts.missing.map(Err),
                                _ => return None,
                            };
                        let Some(allele) = allele else { continue };
                        let buf = &mut buffers[task_i];
                        let dst = &mut buf[patch.idx..patch.idx + patch.rlen];
                        match allele {
                            Ok((bases, case_flags)) => {
                                copy_alt_with_case_flags(dst, bases, patch.to_upper, case_flags);
                                if let Some(mark) = self.opts.mark_snv {
                                    mark_snv_in_place(patch.ref_allele, dst, mark);
                                }
                            }
                            Err(missing) => dst.fill(missing),
                        }
                    }
                }
            }
        }

        let mut remaining = vec![0usize; buffers.len()];
        for &(_, batch_idx) in &output_map {
            remaining[batch_idx] += 1;
        }
        let mut out = Vec::with_capacity(output_map.len());
        for (input_idx, batch_idx) in output_map {
            let src_task = &tasks[input_idx];
            remaining[batch_idx] -= 1;
            let seq = if remaining[batch_idx] == 0 {
                std::mem::take(&mut buffers[batch_idx])
            } else {
                buffers[batch_idx].clone()
            };
            out.push((
                input_idx,
                ConsensusResult {
                    gene_id: src_task.gene_id.clone(),
                    sample: src_task.sample.clone(),
                    haplotype: src_task.haplotype.clone(),
                    seq,
                    chain: None,
                    error: None,
                },
            ));
        }
        Some(out)
    }

    fn try_run_biallelic_phased_alt_only_batch(
        &self,
        tasks: &[ConsensusTask],
        base_ref: &[u8],
        output_map: &[(usize, usize)],
        task_by_hap_sample: &[Vec<usize>; 2],
        active_words_by_hap: &[Vec<u64>; 2],
        active_word_indices_by_hap: &[Vec<usize>; 2],
        patches: &[SameLenBatchPatch<'_>],
    ) -> Option<Vec<(usize, ConsensusResult)>> {
        if self.opts.absent.is_some() || self.opts.missing.is_some() || output_map.is_empty() {
            return None;
        }
        let n_buffers = output_map
            .iter()
            .map(|&(_, batch_idx)| batch_idx)
            .max()
            .map(|idx| idx + 1)?;
        let mut buffers: Vec<Vec<u8>> = (0..n_buffers).map(|_| base_ref.to_vec()).collect();
        for patch in patches {
            match patch {
                SameLenBatchPatch::RefOnly(_) => continue,
                SameLenBatchPatch::Snp1(patch) => {
                    for hap in 1..=2 {
                        let task_by_sample = &task_by_hap_sample[hap - 1];
                        let active_words = &active_words_by_hap[hap - 1];
                        let active_word_indices = &active_word_indices_by_hap[hap - 1];
                        let words = patch.gt_bits.alt_words_for_hap(hap)?;
                        for &word_idx in active_word_indices {
                            let mut bits = *words.get(word_idx)? & active_words[word_idx];
                            while bits != 0 {
                                let bit_idx = bits.trailing_zeros() as usize;
                                let sample_idx = word_idx * 64 + bit_idx;
                                let task_idx = task_by_sample[sample_idx];
                                debug_assert_ne!(task_idx, usize::MAX);
                                if task_idx != usize::MAX {
                                    buffers[task_idx][patch.idx] = patch.alt_out;
                                }
                                bits &= bits - 1;
                            }
                        }
                    }
                }
                SameLenBatchPatch::Mnp(patch) => {
                    for hap in 1..=2 {
                        let task_by_sample = &task_by_hap_sample[hap - 1];
                        let active_words = &active_words_by_hap[hap - 1];
                        let active_word_indices = &active_word_indices_by_hap[hap - 1];
                        let words = patch.gt_bits.alt_words_for_hap(hap)?;
                        for &word_idx in active_word_indices {
                            let mut bits = *words.get(word_idx)? & active_words[word_idx];
                            while bits != 0 {
                                let bit_idx = bits.trailing_zeros() as usize;
                                let sample_idx = word_idx * 64 + bit_idx;
                                let task_idx = task_by_sample[sample_idx];
                                debug_assert_ne!(task_idx, usize::MAX);
                                if task_idx != usize::MAX {
                                    let buf = &mut buffers[task_idx];
                                    let dst = &mut buf[patch.idx..patch.idx + patch.rlen];
                                    copy_alt_with_case_flags(
                                        dst,
                                        patch.alt,
                                        patch.to_upper,
                                        patch.alt_case_flags,
                                    );
                                    if let Some(mark) = self.opts.mark_snv {
                                        mark_snv_in_place(patch.ref_allele, dst, mark);
                                    }
                                }
                                bits &= bits - 1;
                            }
                        }
                    }
                }
            }
        }

        let mut remaining = vec![0usize; buffers.len()];
        for &(_, batch_idx) in output_map {
            remaining[batch_idx] += 1;
        }
        let mut out = Vec::with_capacity(output_map.len());
        for &(input_idx, batch_idx) in output_map {
            let src_task = &tasks[input_idx];
            remaining[batch_idx] -= 1;
            let seq = if remaining[batch_idx] == 0 {
                std::mem::take(&mut buffers[batch_idx])
            } else {
                buffers[batch_idx].clone()
            };
            out.push((
                input_idx,
                ConsensusResult {
                    gene_id: src_task.gene_id.clone(),
                    sample: src_task.sample.clone(),
                    haplotype: src_task.haplotype.clone(),
                    seq,
                    chain: None,
                    error: None,
                },
            ));
        }
        Some(out)
    }

    fn try_run_biallelic_phased_missing_batch(
        &self,
        tasks: &[ConsensusTask],
        base_ref: &[u8],
        output_map: &[(usize, usize)],
        task_by_hap_sample: &[Vec<usize>; 2],
        active_words_by_hap: &[Vec<u64>; 2],
        active_word_indices_by_hap: &[Vec<usize>; 2],
        patches: &[SameLenBatchPatch<'_>],
    ) -> Option<Vec<(usize, ConsensusResult)>> {
        let missing = self.opts.missing?;
        if self.opts.absent.is_some() || output_map.is_empty() {
            return None;
        }
        let n_buffers = output_map
            .iter()
            .map(|&(_, batch_idx)| batch_idx)
            .max()
            .map(|idx| idx + 1)?;
        let mut buffers: Vec<Vec<u8>> = (0..n_buffers).map(|_| base_ref.to_vec()).collect();
        for patch in patches {
            match patch {
                SameLenBatchPatch::RefOnly(_) => continue,
                SameLenBatchPatch::Snp1(patch) => {
                    for hap in 1..=2 {
                        let task_by_sample = &task_by_hap_sample[hap - 1];
                        let active_words = &active_words_by_hap[hap - 1];
                        let active_word_indices = &active_word_indices_by_hap[hap - 1];
                        let alt_words = patch.gt_bits.alt_words_for_hap(hap)?;
                        let missing_words = patch.gt_bits.missing_words();
                        for &word_idx in active_word_indices {
                            let active = active_words[word_idx];
                            let mut alt_bits = *alt_words.get(word_idx)? & active;
                            while alt_bits != 0 {
                                let bit_idx = alt_bits.trailing_zeros() as usize;
                                let sample_idx = word_idx * 64 + bit_idx;
                                let task_idx = task_by_sample[sample_idx];
                                debug_assert_ne!(task_idx, usize::MAX);
                                if task_idx != usize::MAX {
                                    buffers[task_idx][patch.idx] = patch.alt_out;
                                }
                                alt_bits &= alt_bits - 1;
                            }

                            let mut missing_bits = *missing_words.get(word_idx)? & active;
                            while missing_bits != 0 {
                                let bit_idx = missing_bits.trailing_zeros() as usize;
                                let sample_idx = word_idx * 64 + bit_idx;
                                let task_idx = task_by_sample[sample_idx];
                                debug_assert_ne!(task_idx, usize::MAX);
                                if task_idx != usize::MAX {
                                    buffers[task_idx][patch.idx] = missing;
                                }
                                missing_bits &= missing_bits - 1;
                            }
                        }
                    }
                }
                SameLenBatchPatch::Mnp(patch) => {
                    for hap in 1..=2 {
                        let task_by_sample = &task_by_hap_sample[hap - 1];
                        let active_words = &active_words_by_hap[hap - 1];
                        let active_word_indices = &active_word_indices_by_hap[hap - 1];
                        let alt_words = patch.gt_bits.alt_words_for_hap(hap)?;
                        let missing_words = patch.gt_bits.missing_words();
                        for &word_idx in active_word_indices {
                            let active = active_words[word_idx];
                            let mut alt_bits = *alt_words.get(word_idx)? & active;
                            while alt_bits != 0 {
                                let bit_idx = alt_bits.trailing_zeros() as usize;
                                let sample_idx = word_idx * 64 + bit_idx;
                                let task_idx = task_by_sample[sample_idx];
                                debug_assert_ne!(task_idx, usize::MAX);
                                if task_idx != usize::MAX {
                                    let buf = &mut buffers[task_idx];
                                    let dst = &mut buf[patch.idx..patch.idx + patch.rlen];
                                    copy_alt_with_case_flags(
                                        dst,
                                        patch.alt,
                                        patch.to_upper,
                                        patch.alt_case_flags,
                                    );
                                    if let Some(mark) = self.opts.mark_snv {
                                        mark_snv_in_place(patch.ref_allele, dst, mark);
                                    }
                                }
                                alt_bits &= alt_bits - 1;
                            }

                            let mut missing_bits = *missing_words.get(word_idx)? & active;
                            while missing_bits != 0 {
                                let bit_idx = missing_bits.trailing_zeros() as usize;
                                let sample_idx = word_idx * 64 + bit_idx;
                                let task_idx = task_by_sample[sample_idx];
                                debug_assert_ne!(task_idx, usize::MAX);
                                if task_idx != usize::MAX {
                                    let buf = &mut buffers[task_idx];
                                    buf[patch.idx..patch.idx + patch.rlen].fill(missing);
                                }
                                missing_bits &= missing_bits - 1;
                            }
                        }
                    }
                }
            }
        }

        let mut remaining = vec![0usize; buffers.len()];
        for &(_, batch_idx) in output_map {
            remaining[batch_idx] += 1;
        }
        let mut out = Vec::with_capacity(output_map.len());
        for &(input_idx, batch_idx) in output_map {
            let src_task = &tasks[input_idx];
            remaining[batch_idx] -= 1;
            let seq = if remaining[batch_idx] == 0 {
                std::mem::take(&mut buffers[batch_idx])
            } else {
                buffers[batch_idx].clone()
            };
            out.push((
                input_idx,
                ConsensusResult {
                    gene_id: src_task.gene_id.clone(),
                    sample: src_task.sample.clone(),
                    haplotype: src_task.haplotype.clone(),
                    seq,
                    chain: None,
                    error: None,
                },
            ));
        }
        Some(out)
    }

    fn try_run_biallelic_phased_absent_batch(
        &self,
        tasks: &[ConsensusTask],
        base_ref: &[u8],
        output_map: &[(usize, usize)],
        task_by_hap_sample: &[Vec<usize>; 2],
        active_words_by_hap: &[Vec<u64>; 2],
        active_word_indices_by_hap: &[Vec<usize>; 2],
        patches: &[SameLenBatchPatch<'_>],
    ) -> Option<Vec<(usize, ConsensusResult)>> {
        let absent = self.opts.absent?;
        if output_map.is_empty() {
            return None;
        }
        let n_buffers = output_map
            .iter()
            .map(|&(_, batch_idx)| batch_idx)
            .max()
            .map(|idx| idx + 1)?;
        let mut buffers: Vec<Vec<u8>> = (0..n_buffers)
            .map(|_| vec![absent; base_ref.len()])
            .collect();

        for patch in patches {
            match patch {
                SameLenBatchPatch::RefOnly(patch) => {
                    for buf in &mut buffers {
                        if patch.rlen == 1 {
                            buf[patch.idx] = patch.ref_out;
                        } else {
                            let dst = &mut buf[patch.idx..patch.idx + patch.rlen];
                            copy_alt_with_case_flags(
                                dst,
                                patch.ref_allele,
                                patch.to_upper,
                                patch.ref_case_flags,
                            );
                        }
                    }
                }
                SameLenBatchPatch::Snp1(patch) => {
                    for hap in 1..=2 {
                        let task_by_sample = &task_by_hap_sample[hap - 1];
                        let active_words = &active_words_by_hap[hap - 1];
                        let active_word_indices = &active_word_indices_by_hap[hap - 1];
                        let alt_words = patch.gt_bits.alt_words_for_hap(hap)?;
                        let missing_words = patch.gt_bits.missing_words();
                        for &word_idx in active_word_indices {
                            let active = active_words[word_idx];
                            let alt_bits = *alt_words.get(word_idx)? & active;
                            let missing_bits = *missing_words.get(word_idx)? & active;

                            let mut ref_bits = active & !(alt_bits | missing_bits);
                            while ref_bits != 0 {
                                let bit_idx = ref_bits.trailing_zeros() as usize;
                                let sample_idx = word_idx * 64 + bit_idx;
                                let task_idx = task_by_sample[sample_idx];
                                debug_assert_ne!(task_idx, usize::MAX);
                                if task_idx != usize::MAX {
                                    buffers[task_idx][patch.idx] = patch.ref_out;
                                }
                                ref_bits &= ref_bits - 1;
                            }

                            let mut alt_bits = alt_bits;
                            while alt_bits != 0 {
                                let bit_idx = alt_bits.trailing_zeros() as usize;
                                let sample_idx = word_idx * 64 + bit_idx;
                                let task_idx = task_by_sample[sample_idx];
                                debug_assert_ne!(task_idx, usize::MAX);
                                if task_idx != usize::MAX {
                                    buffers[task_idx][patch.idx] = patch.alt_out;
                                }
                                alt_bits &= alt_bits - 1;
                            }

                            if let Some(missing) = self.opts.missing {
                                let mut missing_bits = missing_bits;
                                while missing_bits != 0 {
                                    let bit_idx = missing_bits.trailing_zeros() as usize;
                                    let sample_idx = word_idx * 64 + bit_idx;
                                    let task_idx = task_by_sample[sample_idx];
                                    debug_assert_ne!(task_idx, usize::MAX);
                                    if task_idx != usize::MAX {
                                        buffers[task_idx][patch.idx] = missing;
                                    }
                                    missing_bits &= missing_bits - 1;
                                }
                            }
                        }
                    }
                }
                SameLenBatchPatch::Mnp(patch) => {
                    for hap in 1..=2 {
                        let task_by_sample = &task_by_hap_sample[hap - 1];
                        let active_words = &active_words_by_hap[hap - 1];
                        let active_word_indices = &active_word_indices_by_hap[hap - 1];
                        let alt_words = patch.gt_bits.alt_words_for_hap(hap)?;
                        let missing_words = patch.gt_bits.missing_words();
                        for &word_idx in active_word_indices {
                            let active = active_words[word_idx];
                            let alt_bits = *alt_words.get(word_idx)? & active;
                            let missing_bits = *missing_words.get(word_idx)? & active;

                            let mut ref_bits = active & !(alt_bits | missing_bits);
                            while ref_bits != 0 {
                                let bit_idx = ref_bits.trailing_zeros() as usize;
                                let sample_idx = word_idx * 64 + bit_idx;
                                let task_idx = task_by_sample[sample_idx];
                                debug_assert_ne!(task_idx, usize::MAX);
                                if task_idx != usize::MAX {
                                    let buf = &mut buffers[task_idx];
                                    let dst = &mut buf[patch.idx..patch.idx + patch.rlen];
                                    copy_alt_with_case_flags(
                                        dst,
                                        patch.ref_allele,
                                        patch.to_upper,
                                        patch.ref_case_flags,
                                    );
                                }
                                ref_bits &= ref_bits - 1;
                            }

                            let mut alt_bits = alt_bits;
                            while alt_bits != 0 {
                                let bit_idx = alt_bits.trailing_zeros() as usize;
                                let sample_idx = word_idx * 64 + bit_idx;
                                let task_idx = task_by_sample[sample_idx];
                                debug_assert_ne!(task_idx, usize::MAX);
                                if task_idx != usize::MAX {
                                    let buf = &mut buffers[task_idx];
                                    let dst = &mut buf[patch.idx..patch.idx + patch.rlen];
                                    copy_alt_with_case_flags(
                                        dst,
                                        patch.alt,
                                        patch.to_upper,
                                        patch.alt_case_flags,
                                    );
                                    if let Some(mark) = self.opts.mark_snv {
                                        mark_snv_in_place(patch.ref_allele, dst, mark);
                                    }
                                }
                                alt_bits &= alt_bits - 1;
                            }

                            if let Some(missing) = self.opts.missing {
                                let mut missing_bits = missing_bits;
                                while missing_bits != 0 {
                                    let bit_idx = missing_bits.trailing_zeros() as usize;
                                    let sample_idx = word_idx * 64 + bit_idx;
                                    let task_idx = task_by_sample[sample_idx];
                                    debug_assert_ne!(task_idx, usize::MAX);
                                    if task_idx != usize::MAX {
                                        let buf = &mut buffers[task_idx];
                                        buf[patch.idx..patch.idx + patch.rlen].fill(missing);
                                    }
                                    missing_bits &= missing_bits - 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        let mut remaining = vec![0usize; buffers.len()];
        for &(_, batch_idx) in output_map {
            remaining[batch_idx] += 1;
        }
        let mut out = Vec::with_capacity(output_map.len());
        for &(input_idx, batch_idx) in output_map {
            let src_task = &tasks[input_idx];
            remaining[batch_idx] -= 1;
            let seq = if remaining[batch_idx] == 0 {
                std::mem::take(&mut buffers[batch_idx])
            } else {
                buffers[batch_idx].clone()
            };
            out.push((
                input_idx,
                ConsensusResult {
                    gene_id: src_task.gene_id.clone(),
                    sample: src_task.sample.clone(),
                    haplotype: src_task.haplotype.clone(),
                    seq,
                    chain: None,
                    error: None,
                },
            ));
        }
        Some(out)
    }

    fn plan_options(&self) -> PlanOptions {
        PlanOptions {
            chain: self.opts.chain,
            mark_del: self.opts.mark_del.is_some(),
            mark_ins: self.opts.mark_ins.is_some(),
            mark_snv: self.opts.mark_snv.is_some(),
            mask: self.mask.is_some(),
            mask_skips_variants: self
                .mask
                .as_ref()
                .is_some_and(|mask| mask.with.skips_variants()),
            mask_overlaps_variant: false,
            absent: self.opts.absent.is_some(),
        }
    }

    fn plan_options_for_records(&self, chr: &str, records: &RecordSet<'_>) -> PlanOptions {
        let mut opts = self.plan_options();
        if opts.mask_skips_variants {
            opts.mask_overlaps_variant = self.mask.as_ref().is_some_and(|mask| {
                records
                    .iter()
                    .any(|rec| mask.overlaps(chr, rec.pos, rec.ref_end()))
            });
        }
        opts
    }

    /// Whether chain output is enabled (read by the PyO3 layer).
    pub fn opts_chain(&self) -> bool {
        self.opts.chain
    }

    /// Launch a lazy iterator over results using a producer-consumer model.
    ///
    /// A rayon thread pool (`threads` workers) consumes region groups from an
    /// internal queue and pushes results to a bounded completion queue (capacity =
    /// `prefetch_steps`, giving backpressure). The caller drives consumption via
    /// [`ConsensusIter::next`], which blocks on the completion queue and, upon
    /// receiving a completed group, submits more work to maintain the in-flight
    /// invariant = `prefetch_steps`.
    ///
    /// * `prefetch_steps` — 0 = no prefetch (next submits 1 and blocks on it);
    ///   N>0 = keep N region groups in flight at all times.
    /// * `warmup` — if true, submit the first `prefetch_steps` groups now;
    ///   otherwise submission is deferred to the first `next()` call.
    /// * `ordered` — true = yield results in input order (results buffered by
    ///   index; parallelism capped by the slowest in-order task); false = yield
    ///   in completion order, each result carrying its input `idx`.
    /// * `threads` — rayon pool size for this call (built and torn down per
    ///   iterator; the engine holds no pool of its own).
    pub fn consensus_iter(
        &self,
        tasks: Vec<ConsensusTask>,
        prefetch_steps: usize,
        warmup: bool,
        ordered: bool,
        threads: usize,
    ) -> ConsensusIter {
        let n = tasks.len();
        let nthr = threads.max(1);
        // Completion queue capacity = prefetch window (+1 headroom so a worker
        // can always deliver one result even when next() hasn't drained yet).
        let cap = prefetch_steps.max(1);
        let (res_tx, res_rx): (Sender<GroupResult>, Receiver<GroupResult>) = bounded(cap);
        // Workers pull region-group indices; prefetch logic controls how many
        // groups are outstanding, not the raw task count.
        let (group_tx, group_rx): (Sender<usize>, Receiver<usize>) = bounded(nthr.max(1));
        let groups = Arc::new(group_tasks(&tasks));
        let n_groups = groups.len();
        let tasks = Arc::new(tasks);

        let iter = ConsensusIter {
            tasks: tasks.clone(),
            groups: groups.clone(),
            group_tx: group_tx.clone(),
            res_rx,
            n,
            n_groups,
            submitted_groups: 0,
            completed_groups: 0,
            returned: 0,
            prefetch: prefetch_steps,
            ordered,
            // ordered-mode buffer: idx → result, drained in ascending idx order
            pending: HashMap::new(),
            // unordered-mode buffer: group results already received, drained FIFO
            ready: VecDeque::new(),
            next_yield: 0,
            closed: false,
        };

        if n == 0 {
            // No tasks; workers would immediately exit. Return a drained iter.
            return ConsensusIter {
                closed: true,
                ..iter
            };
        }

        // Spawn the worker pool on a detached thread. It owns `group_rx` and
        // `res_tx`; when `group_tx` is dropped (iterator dropped) the group queue
        // closes and workers exit, which drops `res_tx` and unblocks `next()`.
        let engine = self.clone();
        let pool = thread_pool(nthr);
        std::thread::spawn(move || {
            pool.scope(|scope| {
                for _ in 0..nthr {
                    let group_rx = group_rx.clone();
                    let res_tx = res_tx.clone();
                    let tasks = tasks.clone();
                    let groups = groups.clone();
                    let engine = engine.clone();
                    scope.spawn(move |_| {
                        while let Ok(gidx) = group_rx.recv() {
                            let results = engine.run_group(&tasks, &groups[gidx]);
                            let _ = res_tx.send(GroupResult { results });
                        }
                    });
                }
                // Drop the supervisor's sender clone; the iterator's sender
                // controls lifetime and closes the queue when dropped.
                drop(group_tx);
            });
        });

        let mut iter = iter;
        if warmup {
            iter.submit_up_to();
        }
        iter
    }

    /// Cheap clone for sharing across worker threads (Arc internals).
    fn clone_shallow(&self) -> ConsensusEngine {
        self.clone()
    }
}

/// A completed region group. Each entry is tagged with its original task index.
struct GroupResult {
    results: Vec<(usize, ConsensusResult)>,
}

/// Lazy iterator over consensus results, fed by a background worker pool.
///
/// The prefetch scheduler lives here: `next()` blocks on the completion queue
/// and, after yielding a result, submits more region groups to keep
/// `submitted_groups - completed_groups == prefetch` in flight.
pub struct ConsensusIter {
    /// Held to keep the task data alive for the worker pool's Arc clones.
    #[allow(dead_code)]
    tasks: Arc<Vec<ConsensusTask>>,
    /// Held for the same reason; workers consume group indices and dereference
    /// this shared table.
    #[allow(dead_code)]
    groups: Arc<Vec<TaskGroup>>,
    group_tx: Sender<usize>,
    res_rx: Receiver<GroupResult>,
    n: usize,
    n_groups: usize,
    submitted_groups: usize,
    completed_groups: usize,
    returned: usize,
    prefetch: usize,
    ordered: bool,
    /// ordered mode: results received ahead of their yield index.
    pending: HashMap<usize, ConsensusResult>,
    /// unordered mode: completed group results waiting to be yielded.
    ready: VecDeque<(usize, ConsensusResult)>,
    next_yield: usize,
    closed: bool,
}

impl ConsensusIter {
    /// Submit region groups until the prefetch window is full.
    fn submit_up_to(&mut self) {
        if self.closed {
            return;
        }
        while self.submitted_groups < self.n_groups
            && self.submitted_groups - self.completed_groups < self.prefetch
        {
            let gidx = self.submitted_groups;
            if self.group_tx.send(gidx).is_err() {
                self.closed = true;
                break;
            }
            self.submitted_groups += 1;
        }
        // When prefetch == 0 (no prefetch), submit exactly one and let next()
        // block on it; the next group is only submitted after completion.
        if self.prefetch == 0
            && self.submitted_groups == self.completed_groups
            && self.submitted_groups < self.n_groups
        {
            let gidx = self.submitted_groups;
            if self.group_tx.send(gidx).is_ok() {
                self.submitted_groups += 1;
            } else {
                self.closed = true;
            }
        }
    }

    /// Block until the next result is available, then submit follow-up tasks
    /// to maintain the prefetch invariant. Returns None when exhausted.
    ///
    /// Returns `(idx, result)` where `idx` is the task's input position —
    /// needed for unordered mode so the caller can re-pair 1pIu/2pIu.
    ///
    /// This is GIL-free on the Rust side; the PyO3 wrapper releases the GIL
    /// around this call.
    pub fn next_blocking(&mut self) -> Option<(usize, ConsensusResult)> {
        if self.returned >= self.n {
            return None;
        }

        // Ensure at least one group is in flight before blocking on a result.
        // The first call has submitted nothing yet (no warmup); without this,
        // next() would block on res_rx while no worker has work to run.
        self.submit_up_to();

        let (idx, result) = if self.ordered {
            // Yield strictly in input order; buffer out-of-order arrivals.
            if let Some(r) = self.pending.remove(&self.next_yield) {
                (self.next_yield, r)
            } else {
                // Need to receive until next_yield arrives.
                loop {
                    match self.res_rx.recv() {
                        Ok(group) => {
                            self.completed_groups += 1;
                            self.submit_up_to();
                            let mut next = None;
                            for (got_idx, got_result) in group.results {
                                if got_idx == self.next_yield {
                                    next = Some((got_idx, got_result));
                                } else {
                                    self.pending.insert(got_idx, got_result);
                                }
                            }
                            if let Some(next) = next {
                                break next;
                            }
                        }
                        Err(_) => {
                            self.closed = true;
                            return self
                                .pending
                                .remove(&self.next_yield)
                                .map(|r| (self.next_yield, r));
                        }
                    }
                }
            }
        } else {
            // Unordered: yield whatever finishes first.
            if let Some(ready) = self.ready.pop_front() {
                ready
            } else {
                loop {
                    match self.res_rx.recv() {
                        Ok(group) => {
                            self.completed_groups += 1;
                            self.submit_up_to();
                            self.ready.extend(group.results);
                            if let Some(ready) = self.ready.pop_front() {
                                break ready;
                            }
                        }
                        Err(_) => return None,
                    }
                }
            }
        };

        self.returned += 1;
        if self.ordered {
            self.next_yield += 1;
        }
        // Maintain prefetch invariant: we just freed one slot, submit more.
        self.submit_up_to();
        Some((idx, result))
    }
}

impl Iterator for ConsensusIter {
    type Item = (usize, ConsensusResult);

    fn next(&mut self) -> Option<(usize, ConsensusResult)> {
        self.next_blocking()
    }
}

/// Build the SampleMode for a task from its (sample, haplotype) cli strings.
fn build_sample_mode(
    vcf: &VcfStore,
    sample: &Option<String>,
    haplotype: &Option<String>,
    iupac_codes: bool,
) -> SampleMode {
    match (sample, haplotype) {
        (None, None) if iupac_codes => SampleMode::IupacFromRefAlt,
        (None, None) => SampleMode::ApplyAllAlt,
        (Some(s), None) => {
            let idx = vcf.sample_index(s).unwrap_or(-1);
            SampleMode::IupacAllSamples { samples: vec![idx] }
        }
        (Some(s), Some(h)) => {
            let idx = vcf.sample_index(s).unwrap_or(-1);
            let spec = HaplotypeSpec::parse(h).unwrap_or_default();
            SampleMode::SingleSample { idx, spec }
        }
        (None, Some(h)) => {
            // -H without -s: treat as IupacFromRefAlt-ish; bcftools errors, but
            // we fall back to applying all ALT with the spec ignored.
            let _ = h;
            SampleMode::ApplyAllAlt
        }
    }
}

fn parse_haplotype_index_no_alloc(s: &str) -> Option<u32> {
    let bytes = s.as_bytes();
    let mut i = 0usize;
    let mut hap = 0u32;
    while i < bytes.len() && bytes[i].is_ascii_digit() {
        hap = hap.checked_mul(10)?.checked_add((bytes[i] - b'0') as u32)?;
        i += 1;
    }
    if i == 0 || hap == 0 {
        return None;
    }
    if i == bytes.len() {
        return Some(hap);
    }
    if bytes.len() - i == 3
        && bytes[i].eq_ignore_ascii_case(&b'p')
        && bytes[i + 1].eq_ignore_ascii_case(&b'i')
        && bytes[i + 2].eq_ignore_ascii_case(&b'u')
    {
        return Some(hap);
    }
    None
}

fn task_group_key(task: &ConsensusTask) -> TaskGroupKey {
    TaskGroupKey {
        chr: task.chr.clone(),
        start: task.start,
        end: task.end,
        vcf_key: task.vcf_key.clone(),
    }
}

fn task_exec_key(task: &ConsensusTask) -> TaskExecKey {
    TaskExecKey {
        sample: task.sample.clone(),
        haplotype: task.haplotype.clone(),
    }
}

fn has_duplicate_exec_keys(tasks: &[ConsensusTask], indices: &[usize]) -> bool {
    let mut seen: HashSet<(Option<&str>, Option<&str>)> = HashSet::with_capacity(indices.len());
    for &idx in indices {
        let task = &tasks[idx];
        let key = (task.sample.as_deref(), task.haplotype.as_deref());
        if !seen.insert(key) {
            return true;
        }
    }
    false
}

fn group_tasks(tasks: &[ConsensusTask]) -> Vec<TaskGroup> {
    let mut group_index: HashMap<TaskGroupKey, usize> = HashMap::new();
    let mut groups: Vec<TaskGroup> = Vec::new();
    for (idx, task) in tasks.iter().enumerate() {
        let key = task_group_key(task);
        if let Some(&gidx) = group_index.get(&key) {
            groups[gidx].indices.push(idx);
        } else {
            let gidx = groups.len();
            group_index.insert(key.clone(), gidx);
            groups.push(TaskGroup {
                key,
                indices: vec![idx],
            });
        }
    }
    groups
}

impl From<&ConsensusResult> for CachedOutput {
    fn from(result: &ConsensusResult) -> Self {
        CachedOutput {
            seq: result.seq.clone(),
            chain: result.chain.clone(),
            error: result.error.clone(),
        }
    }
}

fn result_from_cached(task: &ConsensusTask, cached: &CachedOutput) -> ConsensusResult {
    ConsensusResult {
        gene_id: task.gene_id.clone(),
        sample: task.sample.clone(),
        haplotype: task.haplotype.clone(),
        seq: cached.seq.clone(),
        chain: cached.chain.clone(),
        error: cached.error.clone(),
    }
}

fn copy_alt_with_case_flags(dst: &mut [u8], alt: &[u8], to_upper: bool, case_flags: u8) {
    debug_assert_eq!(dst.len(), alt.len());
    if to_upper {
        if case_flags & ALLELE_HAS_ASCII_LOWER == 0 {
            dst.copy_from_slice(alt);
            return;
        }
        if alt.len() == 1 {
            dst[0] = alt[0].to_ascii_uppercase();
            return;
        }
        for (d, &src) in dst.iter_mut().zip(alt) {
            *d = src.to_ascii_uppercase();
        }
    } else {
        if case_flags & ALLELE_HAS_ASCII_UPPER == 0 {
            dst.copy_from_slice(alt);
            return;
        }
        if alt.len() == 1 {
            dst[0] = alt[0].to_ascii_lowercase();
            return;
        }
        for (d, &src) in dst.iter_mut().zip(alt) {
            *d = src.to_ascii_lowercase();
        }
    }
}

#[inline]
fn snp1_alt_with_case_and_mark(
    ref_base: u8,
    alt: u8,
    to_upper: bool,
    case_flags: u8,
    mark_snv: Option<u8>,
) -> u8 {
    let mut out = if to_upper {
        if case_flags & ALLELE_HAS_ASCII_LOWER == 0 {
            alt
        } else {
            alt.to_ascii_uppercase()
        }
    } else if case_flags & ALLELE_HAS_ASCII_UPPER == 0 {
        alt
    } else {
        alt.to_ascii_lowercase()
    };

    if let Some(mark) = mark_snv {
        if !ref_base.eq_ignore_ascii_case(&out) {
            out = if mark == TO_UPPER as u8 {
                out.to_ascii_uppercase()
            } else if mark == TO_LOWER as u8 {
                out.to_ascii_lowercase()
            } else {
                mark
            };
        }
    }
    out
}

fn mark_snv_in_place(ref_allele: &[u8], dst: &mut [u8], mark: u8) {
    let n = ref_allele.len().min(dst.len());
    if mark == TO_UPPER as u8 {
        for i in 0..n {
            if !ref_allele[i].eq_ignore_ascii_case(&dst[i]) {
                dst[i] = dst[i].to_ascii_uppercase();
            }
        }
    } else if mark == TO_LOWER as u8 {
        for i in 0..n {
            if !ref_allele[i].eq_ignore_ascii_case(&dst[i]) {
                dst[i] = dst[i].to_ascii_lowercase();
            }
        }
    } else {
        for i in 0..n {
            if !ref_allele[i].eq_ignore_ascii_case(&dst[i]) {
                dst[i] = mark;
            }
        }
    }
}

fn error_result(task: &ConsensusTask, err: String) -> ConsensusResult {
    ConsensusResult {
        gene_id: task.gene_id.clone(),
        sample: task.sample.clone(),
        haplotype: task.haplotype.clone(),
        seq: Vec::new(),
        chain: None,
        error: Some(err),
    }
}

fn thread_pool(n: usize) -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .expect("failed to build rayon pool")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn setup(name: &str) -> (std::path::PathBuf, std::path::PathBuf) {
        let dir = std::env::temp_dir().join(format!("consensus_rs_engine_test_{}", name));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let ref_fa = dir.join("ref.fa");
        std::fs::write(&ref_fa, format!(">chr1\n{}\n", seq)).unwrap();
        let vcf = dir.join("v.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100>\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t2\t.\tC\tG\t.\t.\t.\n",
        )
        .unwrap();
        (ref_fa, vcf)
    }

    fn write_mask_near(path: &std::path::Path, body: &str) -> std::path::PathBuf {
        let mask = path.parent().unwrap().join("mask.bed");
        std::fs::write(&mask, body).unwrap();
        mask
    }

    fn setup_two_regions(name: &str) -> (std::path::PathBuf, std::path::PathBuf) {
        let dir = std::env::temp_dir().join(format!("consensus_rs_engine_test_{}", name));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let ref_fa = dir.join("ref.fa");
        std::fs::write(&ref_fa, format!(">chr1\n{}\n", seq)).unwrap();
        let vcf = dir.join("v.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100>\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t2\t.\tC\tG\t.\t.\t.\n\
             chr1\t10\t.\tC\tA\t.\t.\t.\n",
        )
        .unwrap();
        (ref_fa, vcf)
    }

    fn setup_phased_batch(name: &str) -> (std::path::PathBuf, std::path::PathBuf) {
        let dir = std::env::temp_dir().join(format!("consensus_rs_engine_test_{}", name));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let ref_fa = dir.join("ref.fa");
        std::fs::write(&ref_fa, ">chr1\nACGTACGT\n").unwrap();
        let vcf = dir.join("v.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100>\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
             chr1\t2\t.\tC\tG\t.\t.\t.\tGT\t0|1\t1|0\n\
             chr1\t3\t.\tG\tT\t.\t.\t.\tGT\t1|0\t0|1\n",
        )
        .unwrap();
        (ref_fa, vcf)
    }

    fn setup_phased_batch_sparse_word(name: &str) -> (std::path::PathBuf, std::path::PathBuf) {
        let dir = std::env::temp_dir().join(format!("consensus_rs_engine_test_{}", name));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let ref_fa = dir.join("ref.fa");
        std::fs::write(&ref_fa, ">chr1\nACGTACGT\n").unwrap();
        let vcf = dir.join("v.vcf");

        let mut header = String::from(
            "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100>\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
        );
        for i in 1..=65 {
            header.push_str(&format!("\tS{}", i));
        }
        header.push('\n');

        let mut row = String::from("chr1\t2\t.\tC\tG\t.\t.\t.\tGT");
        for i in 1..=65 {
            row.push('\t');
            if i == 65 {
                row.push_str("0|1");
            } else {
                row.push_str("0|0");
            }
        }
        row.push('\n');

        std::fs::write(&vcf, format!("{}{}", header, row)).unwrap();
        (ref_fa, vcf)
    }

    fn setup_phased_batch_mnp(name: &str) -> (std::path::PathBuf, std::path::PathBuf) {
        let dir = std::env::temp_dir().join(format!("consensus_rs_engine_test_{}", name));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let ref_fa = dir.join("ref.fa");
        std::fs::write(&ref_fa, ">chr1\nACGTACGT\n").unwrap();
        let vcf = dir.join("v.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100>\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
             chr1\t2\t.\tCG\tTT\t.\t.\t.\tGT\t0|1\t1|0\n",
        )
        .unwrap();
        (ref_fa, vcf)
    }

    fn setup_phased_batch_lowercase(name: &str) -> (std::path::PathBuf, std::path::PathBuf) {
        let dir = std::env::temp_dir().join(format!("consensus_rs_engine_test_{}", name));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let ref_fa = dir.join("ref.fa");
        std::fs::write(&ref_fa, ">chr1\nacgtacgt\n").unwrap();
        let vcf = dir.join("v.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100>\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
             chr1\t2\t.\tc\tG\t.\t.\t.\tGT\t0|1\n\
             chr1\t3\t.\tg\tT\t.\t.\t.\tGT\t1|0\n",
        )
        .unwrap();
        (ref_fa, vcf)
    }

    fn setup_phased_batch_ref_only(name: &str) -> (std::path::PathBuf, std::path::PathBuf) {
        let dir = std::env::temp_dir().join(format!("consensus_rs_engine_test_{}", name));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let ref_fa = dir.join("ref.fa");
        std::fs::write(&ref_fa, ">chr1\nACGTACGT\n").unwrap();
        let vcf = dir.join("v.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100>\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
             chr1\t2\t.\tC\tG\t.\t.\t.\tGT\t0|1\t1|0\n\
             chr1\t5\t.\tA\t.\t.\t.\t.\tGT\t0|0\t0|0\n",
        )
        .unwrap();
        (ref_fa, vcf)
    }

    fn setup_phased_batch_missing(name: &str) -> (std::path::PathBuf, std::path::PathBuf) {
        let dir = std::env::temp_dir().join(format!("consensus_rs_engine_test_{}", name));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let ref_fa = dir.join("ref.fa");
        std::fs::write(&ref_fa, ">chr1\nACGTACGT\n").unwrap();
        let vcf = dir.join("v.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100>\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
             chr1\t2\t.\tC\tG\t.\t.\t.\tGT\t./.\t0|1\n",
        )
        .unwrap();
        (ref_fa, vcf)
    }

    fn setup_phased_batch_mnp_missing(name: &str) -> (std::path::PathBuf, std::path::PathBuf) {
        let dir = std::env::temp_dir().join(format!("consensus_rs_engine_test_{}", name));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let ref_fa = dir.join("ref.fa");
        std::fs::write(&ref_fa, ">chr1\nACGTACGT\n").unwrap();
        let vcf = dir.join("v.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.3\n##contig=<ID=chr1,length=100>\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
             chr1\t2\t.\tCG\tTT\t.\t.\t.\tGT\t./.\t0|1\n",
        )
        .unwrap();
        (ref_fa, vcf)
    }

    #[test]
    fn parses_haplotype_index_without_allocating() {
        assert_eq!(parse_haplotype_index_no_alloc("1"), Some(1));
        assert_eq!(parse_haplotype_index_no_alloc("02"), Some(2));
        assert_eq!(parse_haplotype_index_no_alloc("1pIu"), Some(1));
        assert_eq!(parse_haplotype_index_no_alloc("2PIU"), Some(2));
        assert_eq!(parse_haplotype_index_no_alloc("R"), None);
        assert_eq!(parse_haplotype_index_no_alloc("0"), None);
        assert_eq!(parse_haplotype_index_no_alloc("1x"), None);
    }

    #[test]
    fn groups_tasks_by_region_and_vcf_key() {
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "a".into(),
                gene_id: "G0".into(),
                sample: None,
                haplotype: None,
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 9,
                end: 12,
                vcf_key: "a".into(),
                gene_id: "G1".into(),
                sample: None,
                haplotype: None,
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "a".into(),
                gene_id: "G2".into(),
                sample: Some("S1".into()),
                haplotype: Some("1".into()),
            },
        ];
        let groups = group_tasks(&tasks);
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].indices, vec![0, 2]);
        assert_eq!(groups[1].indices, vec![1]);
        assert!(!has_duplicate_exec_keys(&tasks, &groups[0].indices));

        let duplicate = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "a".into(),
                gene_id: "D0".into(),
                sample: Some("S1".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "a".into(),
                gene_id: "D1".into(),
                sample: Some("S1".into()),
                haplotype: Some("1".into()),
            },
        ];
        let duplicate_groups = group_tasks(&duplicate);
        assert!(has_duplicate_exec_keys(
            &duplicate,
            &duplicate_groups[0].indices
        ));
    }

    #[test]
    fn engine_single_task_snp() {
        let (ref_fa, vcf) = setup("single");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let opts = EngineOptions {
            ..Default::default()
        };
        let engine = ConsensusEngine::load(ref_fa, vcf_map, opts).unwrap();
        let task = ConsensusTask {
            chr: "chr1".into(),
            start: 1,
            end: 8,
            vcf_key: "chr1".into(),
            gene_id: "G1".into(),
            sample: None,
            haplotype: None,
        };
        let results = engine.consensus_many(vec![task], 2);
        assert_eq!(results.len(), 1);
        assert!(results[0].error.is_none(), "{:?}", results[0].error);
        // chr1:2 C>G over ACGTACGT -> AGGTACGT
        assert_eq!(results[0].seq, b"AGGTACGT");
        assert_eq!(results[0].gene_id, "G1");
    }

    #[test]
    fn engine_many_tasks_parallel() {
        let (ref_fa, vcf) = setup("many");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let opts = EngineOptions {
            ..Default::default()
        };
        let engine = ConsensusEngine::load(ref_fa, vcf_map, opts).unwrap();
        let tasks: Vec<ConsensusTask> = (0..20)
            .map(|i| ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: format!("G{}", i),
                sample: None,
                haplotype: None,
            })
            .collect();
        let results = engine.consensus_many(tasks, 4);
        assert_eq!(results.len(), 20);
        for r in &results {
            assert!(r.error.is_none(), "{:?}", r.error);
            assert_eq!(r.seq, b"AGGTACGT");
        }
        // results preserve input order
        for (i, r) in results.iter().enumerate() {
            assert_eq!(r.gene_id, format!("G{}", i));
        }
    }

    #[test]
    fn engine_many_groups_preserves_input_order() {
        let (ref_fa, vcf) = setup_two_regions("many_groups");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(ref_fa, vcf_map, EngineOptions::default()).unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "G0".into(),
                sample: None,
                haplotype: None,
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 9,
                end: 12,
                vcf_key: "chr1".into(),
                gene_id: "G1".into(),
                sample: None,
                haplotype: None,
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "G2".into(),
                sample: None,
                haplotype: None,
            },
        ];
        let results = engine.consensus_many(tasks, 4);
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].gene_id, "G0");
        assert_eq!(results[1].gene_id, "G1");
        assert_eq!(results[2].gene_id, "G2");
        assert_eq!(results[0].seq, b"AGGTACGT");
        assert_eq!(results[1].seq, b"AAGT");
        assert_eq!(results[2].seq, b"AGGTACGT");
    }

    #[test]
    fn empty_region_group_fastpath_handles_absent_and_chain() {
        let (ref_fa, vcf) = setup("empty_group_absent_chain");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(
            ref_fa,
            vcf_map,
            EngineOptions {
                absent: Some(b'N'),
                chain: true,
                ..Default::default()
            },
        )
        .unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 20,
                end: 23,
                vcf_key: "chr1".into(),
                gene_id: "E0".into(),
                sample: None,
                haplotype: None,
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 20,
                end: 23,
                vcf_key: "chr1".into(),
                gene_id: "E1".into(),
                sample: None,
                haplotype: None,
            },
        ];

        let results = engine.consensus_many(tasks, 2);

        assert_eq!(results.len(), 2);
        assert_eq!(results[0].seq, b"NNNN");
        assert_eq!(results[1].seq, b"NNNN");
        assert_eq!(
            results[0].chain.as_deref(),
            Some("chain 4 chr1 23 + 19 23 chr1 23 + 19 23 1\n4\n\n")
        );
        assert_eq!(results[1].chain, results[0].chain);
    }

    #[test]
    fn engine_uses_compiled_mask_and_skips_overlapping_variant() {
        let (ref_fa, vcf) = setup("compiled_mask_engine");
        let mask = vcf.parent().unwrap().join("mask.bed");
        std::fs::write(&mask, "chr1\t1\t2\n").unwrap();
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(
            ref_fa,
            vcf_map,
            EngineOptions {
                mask: Some(mask),
                mask_with: crate::mask::MaskWith::Char(b'N'),
                ..Default::default()
            },
        )
        .unwrap();
        let task = ConsensusTask {
            chr: "chr1".into(),
            start: 1,
            end: 8,
            vcf_key: "chr1".into(),
            gene_id: "M0".into(),
            sample: None,
            haplotype: None,
        };

        let results = engine.consensus_many(vec![task], 1);

        assert_eq!(results.len(), 1);
        assert!(results[0].error.is_none(), "{:?}", results[0].error);
        assert_eq!(results[0].seq, b"ANGTACGT");
    }

    #[test]
    fn consensus_iter_ordered_uses_region_groups() {
        let (ref_fa, vcf) = setup_two_regions("iter_ordered_groups");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(ref_fa, vcf_map, EngineOptions::default()).unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "G0".into(),
                sample: None,
                haplotype: None,
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 9,
                end: 12,
                vcf_key: "chr1".into(),
                gene_id: "G1".into(),
                sample: None,
                haplotype: None,
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "G2".into(),
                sample: None,
                haplotype: None,
            },
        ];

        let results: Vec<_> = engine.consensus_iter(tasks, 1, true, true, 2).collect();
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].0, 0);
        assert_eq!(results[1].0, 1);
        assert_eq!(results[2].0, 2);
        assert_eq!(results[0].1.seq, b"AGGTACGT");
        assert_eq!(results[1].1.seq, b"AAGT");
        assert_eq!(results[2].1.seq, b"AGGTACGT");
    }

    #[test]
    fn consensus_iter_unordered_prefetch_zero_drains_group_results() {
        let (ref_fa, vcf) = setup_two_regions("iter_unordered_prefetch_zero");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(ref_fa, vcf_map, EngineOptions::default()).unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "G0".into(),
                sample: None,
                haplotype: None,
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 9,
                end: 12,
                vcf_key: "chr1".into(),
                gene_id: "G1".into(),
                sample: None,
                haplotype: None,
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "G2".into(),
                sample: None,
                haplotype: None,
            },
        ];

        let mut results: Vec<_> = engine.consensus_iter(tasks, 0, false, false, 2).collect();
        results.sort_by_key(|(idx, _)| *idx);
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].0, 0);
        assert_eq!(results[1].0, 1);
        assert_eq!(results[2].0, 2);
        assert_eq!(results[0].1.seq, b"AGGTACGT");
        assert_eq!(results[1].1.seq, b"AAGT");
        assert_eq!(results[2].1.seq, b"AGGTACGT");
    }

    #[test]
    fn biallelic_phased_batch_lane_patches_group_outputs() {
        let (ref_fa, vcf) = setup_phased_batch("phased_batch");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(ref_fa, vcf_map, EngineOptions::default()).unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H1".into(),
                sample: Some("S1".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H2".into(),
                sample: Some("S1".into()),
                haplotype: Some("2".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H1".into(),
                sample: Some("S2".into()),
                haplotype: Some("1pIu".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H2".into(),
                sample: Some("S2".into()),
                haplotype: Some("2pIu".into()),
            },
        ];

        let groups = group_tasks(&tasks);
        let vcf = engine.vcfs.get("chr1").unwrap();
        let ref_seq = engine.ref_index.fetch_1based("chr1", 1, 8).unwrap();
        let (records, plan) = vcf.plan_query_set(
            "chr1",
            0,
            7,
            engine.opts.regions_overlap,
            engine.plan_options(),
        );
        let batch = engine
            .try_run_biallelic_phased_batch(&tasks, &groups[0], vcf, &ref_seq, 0, &records, &plan)
            .expect("batch lane should accept biallelic phased same-len group");
        let mut by_idx = vec![Vec::new(); batch.len()];
        for (idx, result) in batch {
            by_idx[idx] = result.seq;
        }
        assert_eq!(by_idx[0], b"ACTTACGT");
        assert_eq!(by_idx[1], b"AGGTACGT");
        assert_eq!(by_idx[2], b"AGGTACGT");
        assert_eq!(by_idx[3], b"ACTTACGT");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"ACTTACGT");
        assert_eq!(results[1].seq, b"AGGTACGT");
        assert_eq!(results[2].seq, b"AGGTACGT");
        assert_eq!(results[3].seq, b"ACTTACGT");
    }

    #[test]
    fn biallelic_phased_batch_lane_handles_sparse_active_words() {
        let (ref_fa, vcf) = setup_phased_batch_sparse_word("phased_batch_sparse_word");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(ref_fa, vcf_map, EngineOptions::default()).unwrap();
        let tasks = vec![ConsensusTask {
            chr: "chr1".into(),
            start: 1,
            end: 8,
            vcf_key: "chr1".into(),
            gene_id: "S65H2".into(),
            sample: Some("S65".into()),
            haplotype: Some("2".into()),
        }];

        let groups = group_tasks(&tasks);
        let vcf = engine.vcfs.get("chr1").unwrap();
        let ref_seq = engine.ref_index.fetch_1based("chr1", 1, 8).unwrap();
        let (records, plan) = vcf.plan_query_set(
            "chr1",
            0,
            7,
            engine.opts.regions_overlap,
            engine.plan_options(),
        );
        let batch = engine
            .try_run_biallelic_phased_batch(&tasks, &groups[0], vcf, &ref_seq, 0, &records, &plan)
            .expect("batch lane should patch requested samples in sparse active words");
        assert_eq!(batch.len(), 1);
        assert_eq!(batch[0].0, 0);
        assert_eq!(batch[0].1.seq, b"AGGTACGT");

        let results = engine.consensus_many(tasks, 1);
        assert_eq!(results[0].seq, b"AGGTACGT");
    }

    #[test]
    fn biallelic_phased_batch_lane_handles_mnp_patch_plan() {
        let (ref_fa, vcf) = setup_phased_batch_mnp("phased_batch_mnp_patch_plan");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(ref_fa, vcf_map, EngineOptions::default()).unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H1".into(),
                sample: Some("S1".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H2".into(),
                sample: Some("S1".into()),
                haplotype: Some("2".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H1".into(),
                sample: Some("S2".into()),
                haplotype: Some("1".into()),
            },
        ];

        let groups = group_tasks(&tasks);
        let vcf = engine.vcfs.get("chr1").unwrap();
        let ref_seq = engine.ref_index.fetch_1based("chr1", 1, 8).unwrap();
        let (records, plan) = vcf.plan_query_set(
            "chr1",
            0,
            7,
            engine.opts.regions_overlap,
            engine.plan_options(),
        );
        let batch = engine
            .try_run_biallelic_phased_batch(&tasks, &groups[0], vcf, &ref_seq, 0, &records, &plan)
            .expect("batch lane should accept biallelic phased MNP records");
        let mut by_idx = vec![Vec::new(); batch.len()];
        for (idx, result) in batch {
            by_idx[idx] = result.seq;
        }
        assert_eq!(by_idx[0], b"ACGTACGT");
        assert_eq!(by_idx[1], b"ATTTACGT");
        assert_eq!(by_idx[2], b"ATTTACGT");

        let results = engine.consensus_many(tasks, 1);
        assert_eq!(results[0].seq, b"ACGTACGT");
        assert_eq!(results[1].seq, b"ATTTACGT");
        assert_eq!(results[2].seq, b"ATTTACGT");
    }

    #[test]
    fn biallelic_phased_batch_lane_handles_duplicate_exec_keys() {
        let (ref_fa, vcf) = setup_phased_batch("phased_batch_duplicate_exec");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(ref_fa, vcf_map, EngineOptions::default()).unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H2_A".into(),
                sample: Some("S1".into()),
                haplotype: Some("2".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H2_B".into(),
                sample: Some("S1".into()),
                haplotype: Some("2".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H1".into(),
                sample: Some("S2".into()),
                haplotype: Some("1".into()),
            },
        ];

        let groups = group_tasks(&tasks);
        let vcf = engine.vcfs.get("chr1").unwrap();
        let ref_seq = engine.ref_index.fetch_1based("chr1", 1, 8).unwrap();
        let records = vcf.query_set("chr1", 0, 7, engine.opts.regions_overlap);
        let plan = plan_region_set(&records, engine.plan_options_for_records("chr1", &records));
        let batch = engine
            .try_run_biallelic_phased_batch(&tasks, &groups[0], vcf, &ref_seq, 0, &records, &plan)
            .expect("batch lane should deduplicate identical sample/hap outputs");
        let mut by_idx = vec![Vec::new(); batch.len()];
        for (idx, result) in batch {
            by_idx[idx] = result.seq;
        }
        assert_eq!(by_idx[0], b"AGGTACGT");
        assert_eq!(by_idx[1], b"AGGTACGT");
        assert_eq!(by_idx[2], b"AGGTACGT");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"AGGTACGT");
        assert_eq!(results[1].seq, b"AGGTACGT");
        assert_eq!(results[2].seq, b"AGGTACGT");
    }

    #[test]
    fn biallelic_phased_batch_lane_accepts_nonoverlapping_char_mask() {
        let (ref_fa, vcf) = setup_phased_batch("phased_batch_char_mask");
        let mask = write_mask_near(&vcf, "chr1\t4\t5\n");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(
            ref_fa,
            vcf_map,
            EngineOptions {
                mask: Some(mask),
                mask_with: crate::mask::MaskWith::Char(b'N'),
                ..Default::default()
            },
        )
        .unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H1".into(),
                sample: Some("S1".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H2".into(),
                sample: Some("S1".into()),
                haplotype: Some("2".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H1".into(),
                sample: Some("S2".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H2".into(),
                sample: Some("S2".into()),
                haplotype: Some("2".into()),
            },
        ];

        let groups = group_tasks(&tasks);
        let vcf = engine.vcfs.get("chr1").unwrap();
        let ref_seq = engine.ref_index.fetch_1based("chr1", 1, 8).unwrap();
        let records = vcf.query_set("chr1", 0, 7, engine.opts.regions_overlap);
        let plan = plan_region_set(&records, engine.plan_options_for_records("chr1", &records));
        let batch = engine
            .try_run_biallelic_phased_batch(&tasks, &groups[0], vcf, &ref_seq, 0, &records, &plan)
            .expect("batch lane should accept char mask that does not overlap variants");
        let mut by_idx = vec![Vec::new(); batch.len()];
        for (idx, result) in batch {
            by_idx[idx] = result.seq;
        }
        assert_eq!(by_idx[0], b"ACTTNCGT");
        assert_eq!(by_idx[1], b"AGGTNCGT");
        assert_eq!(by_idx[2], b"AGGTNCGT");
        assert_eq!(by_idx[3], b"ACTTNCGT");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"ACTTNCGT");
        assert_eq!(results[1].seq, b"AGGTNCGT");
        assert_eq!(results[2].seq, b"AGGTNCGT");
        assert_eq!(results[3].seq, b"ACTTNCGT");
    }

    #[test]
    fn biallelic_phased_batch_lane_accepts_case_mask_overlap() {
        let (ref_fa, vcf) = setup_phased_batch("phased_batch_lc_mask");
        let mask = write_mask_near(&vcf, "chr1\t1\t2\n");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(
            ref_fa,
            vcf_map,
            EngineOptions {
                mask: Some(mask),
                mask_with: crate::mask::MaskWith::Lc,
                ..Default::default()
            },
        )
        .unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H1".into(),
                sample: Some("S1".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H2".into(),
                sample: Some("S1".into()),
                haplotype: Some("2".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H1".into(),
                sample: Some("S2".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H2".into(),
                sample: Some("S2".into()),
                haplotype: Some("2".into()),
            },
        ];

        let groups = group_tasks(&tasks);
        let vcf = engine.vcfs.get("chr1").unwrap();
        let ref_seq = engine.ref_index.fetch_1based("chr1", 1, 8).unwrap();
        let records = vcf.query_set("chr1", 0, 7, engine.opts.regions_overlap);
        let plan = plan_region_set(&records, engine.plan_options_for_records("chr1", &records));
        let batch = engine
            .try_run_biallelic_phased_batch(&tasks, &groups[0], vcf, &ref_seq, 0, &records, &plan)
            .expect("batch lane should accept case masks overlapping variants");
        let mut by_idx = vec![Vec::new(); batch.len()];
        for (idx, result) in batch {
            by_idx[idx] = result.seq;
        }
        assert_eq!(by_idx[0], b"AcTTACGT");
        assert_eq!(by_idx[1], b"AgGTACGT");
        assert_eq!(by_idx[2], b"AgGTACGT");
        assert_eq!(by_idx[3], b"AcTTACGT");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"AcTTACGT");
        assert_eq!(results[1].seq, b"AgGTACGT");
        assert_eq!(results[2].seq, b"AgGTACGT");
        assert_eq!(results[3].seq, b"AcTTACGT");
    }

    #[test]
    fn biallelic_phased_batch_lane_handles_mark_snv() {
        let (ref_fa, vcf) = setup_phased_batch("phased_batch_mark_snv");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(
            ref_fa,
            vcf_map,
            EngineOptions {
                mark_snv: Some(b'X'),
                ..Default::default()
            },
        )
        .unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H1".into(),
                sample: Some("S1".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H2".into(),
                sample: Some("S1".into()),
                haplotype: Some("2".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H1".into(),
                sample: Some("S2".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H2".into(),
                sample: Some("S2".into()),
                haplotype: Some("2".into()),
            },
        ];

        let groups = group_tasks(&tasks);
        let vcf = engine.vcfs.get("chr1").unwrap();
        let ref_seq = engine.ref_index.fetch_1based("chr1", 1, 8).unwrap();
        let (records, plan) = vcf.plan_query_set(
            "chr1",
            0,
            7,
            engine.opts.regions_overlap,
            engine.plan_options(),
        );
        let batch = engine
            .try_run_biallelic_phased_batch(&tasks, &groups[0], vcf, &ref_seq, 0, &records, &plan)
            .expect("batch lane should accept same-len marked SNPs");
        let mut by_idx = vec![Vec::new(); batch.len()];
        for (idx, result) in batch {
            by_idx[idx] = result.seq;
        }
        assert_eq!(by_idx[0], b"ACXTACGT");
        assert_eq!(by_idx[1], b"AXGTACGT");
        assert_eq!(by_idx[2], b"AXGTACGT");
        assert_eq!(by_idx[3], b"ACXTACGT");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"ACXTACGT");
        assert_eq!(results[1].seq, b"AXGTACGT");
        assert_eq!(results[2].seq, b"AXGTACGT");
        assert_eq!(results[3].seq, b"ACXTACGT");
    }

    #[test]
    fn biallelic_phased_batch_lane_syncs_lowercase_ref_case() {
        let (ref_fa, vcf) = setup_phased_batch_lowercase("phased_batch_lowercase");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(ref_fa, vcf_map, EngineOptions::default()).unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H1".into(),
                sample: Some("S1".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H2".into(),
                sample: Some("S1".into()),
                haplotype: Some("2".into()),
            },
        ];

        let groups = group_tasks(&tasks);
        let vcf = engine.vcfs.get("chr1").unwrap();
        let ref_seq = engine.ref_index.fetch_1based("chr1", 1, 8).unwrap();
        let (records, plan) = vcf.plan_query_set(
            "chr1",
            0,
            7,
            engine.opts.regions_overlap,
            engine.plan_options(),
        );
        let batch = engine
            .try_run_biallelic_phased_batch(&tasks, &groups[0], vcf, &ref_seq, 0, &records, &plan)
            .expect("batch lane should accept lowercase same-len phased records");
        let mut by_idx = vec![Vec::new(); batch.len()];
        for (idx, result) in batch {
            by_idx[idx] = result.seq;
        }
        assert_eq!(by_idx[0], b"acttacgt");
        assert_eq!(by_idx[1], b"aggtacgt");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"acttacgt");
        assert_eq!(results[1].seq, b"aggtacgt");
    }

    #[test]
    fn biallelic_phased_batch_lane_handles_absent_fill() {
        let (ref_fa, vcf) = setup_phased_batch("phased_batch_absent");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(
            ref_fa,
            vcf_map,
            EngineOptions {
                absent: Some(b'N'),
                ..Default::default()
            },
        )
        .unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H1".into(),
                sample: Some("S1".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H2".into(),
                sample: Some("S1".into()),
                haplotype: Some("2".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H1".into(),
                sample: Some("S2".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H2".into(),
                sample: Some("S2".into()),
                haplotype: Some("2".into()),
            },
        ];

        let groups = group_tasks(&tasks);
        let vcf = engine.vcfs.get("chr1").unwrap();
        let ref_seq = engine.ref_index.fetch_1based("chr1", 1, 8).unwrap();
        let (records, plan) = vcf.plan_query_set(
            "chr1",
            0,
            7,
            engine.opts.regions_overlap,
            engine.plan_options(),
        );
        let batch = engine
            .try_run_biallelic_phased_batch(&tasks, &groups[0], vcf, &ref_seq, 0, &records, &plan)
            .expect("batch lane should accept absent fill for same-len phased records");
        let mut by_idx = vec![Vec::new(); batch.len()];
        for (idx, result) in batch {
            by_idx[idx] = result.seq;
        }
        assert_eq!(by_idx[0], b"NCTNNNNN");
        assert_eq!(by_idx[1], b"NGGNNNNN");
        assert_eq!(by_idx[2], b"NGGNNNNN");
        assert_eq!(by_idx[3], b"NCTNNNNN");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"NCTNNNNN");
        assert_eq!(results[1].seq, b"NGGNNNNN");
        assert_eq!(results[2].seq, b"NGGNNNNN");
        assert_eq!(results[3].seq, b"NCTNNNNN");
    }

    #[test]
    fn biallelic_phased_batch_lane_handles_absent_ref_only_records() {
        let (ref_fa, vcf) = setup_phased_batch_ref_only("phased_batch_absent_ref_only");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(
            ref_fa,
            vcf_map,
            EngineOptions {
                absent: Some(b'N'),
                ..Default::default()
            },
        )
        .unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H1".into(),
                sample: Some("S1".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H2".into(),
                sample: Some("S1".into()),
                haplotype: Some("2".into()),
            },
        ];

        let groups = group_tasks(&tasks);
        let vcf = engine.vcfs.get("chr1").unwrap();
        let ref_seq = engine.ref_index.fetch_1based("chr1", 1, 8).unwrap();
        let (records, plan) = vcf.plan_query_set(
            "chr1",
            0,
            7,
            engine.opts.regions_overlap,
            engine.plan_options(),
        );
        let batch = engine
            .try_run_biallelic_phased_batch(&tasks, &groups[0], vcf, &ref_seq, 0, &records, &plan)
            .expect("batch lane should write ref-only spans under absent fill");
        let mut by_idx = vec![Vec::new(); batch.len()];
        for (idx, result) in batch {
            by_idx[idx] = result.seq;
        }
        assert_eq!(by_idx[0], b"NCNNANNN");
        assert_eq!(by_idx[1], b"NGNNANNN");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"NCNNANNN");
        assert_eq!(results[1].seq, b"NGNNANNN");
    }

    #[test]
    fn biallelic_phased_batch_lane_handles_mnp_absent_fill() {
        let (ref_fa, vcf) = setup_phased_batch_mnp("phased_batch_mnp_absent");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(
            ref_fa,
            vcf_map,
            EngineOptions {
                absent: Some(b'N'),
                ..Default::default()
            },
        )
        .unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H1".into(),
                sample: Some("S1".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H2".into(),
                sample: Some("S1".into()),
                haplotype: Some("2".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H1".into(),
                sample: Some("S2".into()),
                haplotype: Some("1".into()),
            },
        ];

        let groups = group_tasks(&tasks);
        let vcf = engine.vcfs.get("chr1").unwrap();
        let ref_seq = engine.ref_index.fetch_1based("chr1", 1, 8).unwrap();
        let (records, plan) = vcf.plan_query_set(
            "chr1",
            0,
            7,
            engine.opts.regions_overlap,
            engine.plan_options(),
        );
        let batch = engine
            .try_run_biallelic_phased_batch(&tasks, &groups[0], vcf, &ref_seq, 0, &records, &plan)
            .expect("batch lane should accept absent fill for phased MNP records");
        let mut by_idx = vec![Vec::new(); batch.len()];
        for (idx, result) in batch {
            by_idx[idx] = result.seq;
        }
        assert_eq!(by_idx[0], b"NCGNNNNN");
        assert_eq!(by_idx[1], b"NTTNNNNN");
        assert_eq!(by_idx[2], b"NTTNNNNN");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"NCGNNNNN");
        assert_eq!(results[1].seq, b"NTTNNNNN");
        assert_eq!(results[2].seq, b"NTTNNNNN");
    }

    #[test]
    fn biallelic_phased_batch_lane_handles_missing_char() {
        let (ref_fa, vcf) = setup_phased_batch_missing("phased_batch_missing");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(
            ref_fa,
            vcf_map,
            EngineOptions {
                missing: Some(b'?'),
                ..Default::default()
            },
        )
        .unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H1".into(),
                sample: Some("S1".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H2".into(),
                sample: Some("S1".into()),
                haplotype: Some("2".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H1".into(),
                sample: Some("S2".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H2".into(),
                sample: Some("S2".into()),
                haplotype: Some("2".into()),
            },
        ];

        let groups = group_tasks(&tasks);
        let vcf = engine.vcfs.get("chr1").unwrap();
        let ref_seq = engine.ref_index.fetch_1based("chr1", 1, 8).unwrap();
        let (records, plan) = vcf.plan_query_set(
            "chr1",
            0,
            7,
            engine.opts.regions_overlap,
            engine.plan_options(),
        );
        let batch = engine
            .try_run_biallelic_phased_batch(&tasks, &groups[0], vcf, &ref_seq, 0, &records, &plan)
            .expect("batch lane should accept missing GT with -M");
        let mut by_idx = vec![Vec::new(); batch.len()];
        for (idx, result) in batch {
            by_idx[idx] = result.seq;
        }
        assert_eq!(by_idx[0], b"A?GTACGT");
        assert_eq!(by_idx[1], b"A?GTACGT");
        assert_eq!(by_idx[2], b"ACGTACGT");
        assert_eq!(by_idx[3], b"AGGTACGT");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"A?GTACGT");
        assert_eq!(results[1].seq, b"A?GTACGT");
        assert_eq!(results[2].seq, b"ACGTACGT");
        assert_eq!(results[3].seq, b"AGGTACGT");
    }

    #[test]
    fn biallelic_phased_batch_lane_handles_mnp_missing_char() {
        let (ref_fa, vcf) = setup_phased_batch_mnp_missing("phased_batch_mnp_missing");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(
            ref_fa,
            vcf_map,
            EngineOptions {
                missing: Some(b'?'),
                ..Default::default()
            },
        )
        .unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H1".into(),
                sample: Some("S1".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H2".into(),
                sample: Some("S1".into()),
                haplotype: Some("2".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H1".into(),
                sample: Some("S2".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H2".into(),
                sample: Some("S2".into()),
                haplotype: Some("2".into()),
            },
        ];

        let groups = group_tasks(&tasks);
        let vcf = engine.vcfs.get("chr1").unwrap();
        let ref_seq = engine.ref_index.fetch_1based("chr1", 1, 8).unwrap();
        let (records, plan) = vcf.plan_query_set(
            "chr1",
            0,
            7,
            engine.opts.regions_overlap,
            engine.plan_options(),
        );
        let batch = engine
            .try_run_biallelic_phased_batch(&tasks, &groups[0], vcf, &ref_seq, 0, &records, &plan)
            .expect("batch lane should accept MNP missing GT with -M");
        let mut by_idx = vec![Vec::new(); batch.len()];
        for (idx, result) in batch {
            by_idx[idx] = result.seq;
        }
        assert_eq!(by_idx[0], b"A??TACGT");
        assert_eq!(by_idx[1], b"A??TACGT");
        assert_eq!(by_idx[2], b"ACGTACGT");
        assert_eq!(by_idx[3], b"ATTTACGT");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"A??TACGT");
        assert_eq!(results[1].seq, b"A??TACGT");
        assert_eq!(results[2].seq, b"ACGTACGT");
        assert_eq!(results[3].seq, b"ATTTACGT");
    }

    #[test]
    fn biallelic_phased_batch_lane_handles_mnp_missing_char_with_absent_fill() {
        let (ref_fa, vcf) = setup_phased_batch_mnp_missing("phased_batch_mnp_missing_absent");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(
            ref_fa,
            vcf_map,
            EngineOptions {
                missing: Some(b'?'),
                absent: Some(b'N'),
                ..Default::default()
            },
        )
        .unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H1".into(),
                sample: Some("S1".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H1".into(),
                sample: Some("S2".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H2".into(),
                sample: Some("S2".into()),
                haplotype: Some("2".into()),
            },
        ];

        let groups = group_tasks(&tasks);
        let vcf = engine.vcfs.get("chr1").unwrap();
        let ref_seq = engine.ref_index.fetch_1based("chr1", 1, 8).unwrap();
        let (records, plan) = vcf.plan_query_set(
            "chr1",
            0,
            7,
            engine.opts.regions_overlap,
            engine.plan_options(),
        );
        let batch = engine
            .try_run_biallelic_phased_batch(&tasks, &groups[0], vcf, &ref_seq, 0, &records, &plan)
            .expect("batch lane should accept MNP -M with absent fill");
        let mut by_idx = vec![Vec::new(); batch.len()];
        for (idx, result) in batch {
            by_idx[idx] = result.seq;
        }
        assert_eq!(by_idx[0], b"N??NNNNN");
        assert_eq!(by_idx[1], b"NCGNNNNN");
        assert_eq!(by_idx[2], b"NTTNNNNN");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"N??NNNNN");
        assert_eq!(results[1].seq, b"NCGNNNNN");
        assert_eq!(results[2].seq, b"NTTNNNNN");
    }

    #[test]
    fn biallelic_phased_batch_lane_handles_missing_char_with_absent_fill() {
        let (ref_fa, vcf) = setup_phased_batch_missing("phased_batch_missing_absent");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(
            ref_fa,
            vcf_map,
            EngineOptions {
                missing: Some(b'?'),
                absent: Some(b'N'),
                ..Default::default()
            },
        )
        .unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H1".into(),
                sample: Some("S1".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H1".into(),
                sample: Some("S2".into()),
                haplotype: Some("1".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S2H2".into(),
                sample: Some("S2".into()),
                haplotype: Some("2".into()),
            },
        ];

        let groups = group_tasks(&tasks);
        let vcf = engine.vcfs.get("chr1").unwrap();
        let ref_seq = engine.ref_index.fetch_1based("chr1", 1, 8).unwrap();
        let (records, plan) = vcf.plan_query_set(
            "chr1",
            0,
            7,
            engine.opts.regions_overlap,
            engine.plan_options(),
        );
        let batch = engine
            .try_run_biallelic_phased_batch(&tasks, &groups[0], vcf, &ref_seq, 0, &records, &plan)
            .expect("batch lane should accept -M with absent fill");
        let mut by_idx = vec![Vec::new(); batch.len()];
        for (idx, result) in batch {
            by_idx[idx] = result.seq;
        }
        assert_eq!(by_idx[0], b"N?NNNNNN");
        assert_eq!(by_idx[1], b"NCNNNNNN");
        assert_eq!(by_idx[2], b"NGNNNNNN");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"N?NNNNNN");
        assert_eq!(results[1].seq, b"NCNNNNNN");
        assert_eq!(results[2].seq, b"NGNNNNNN");
    }
}
