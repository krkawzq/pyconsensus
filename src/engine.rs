//! engine — ConsensusEngine: holds preprocessed material (ref + vcfs), runs
//! multi-threaded consensus production.
//!
//! (docs/design.md §5.4 / §5.6) Each (region, sample, haplotype) task is
//! independent → natural parallelism. `consensus_many` and `consensus_iter`
//! group identical regions so ref fetch, VCF query, and region planning are
//! amortized across sample/haplotype tasks.
//!
//! This module is PyO3-free; `py.rs` wraps it under the `python` feature.

use crate::apply::{apply_region_planned, ApplyOptions, TO_LOWER, TO_UPPER};
use crate::chain::Chain;
use crate::compiled::RecordFlags;
use crate::haplotype::{HaplotypeSpec, SampleMode};
use crate::planner::RegionPlan;
use crate::ref_index::RefIndex;
use crate::stats::FastPathLane;
use crate::vcf_store::{LoadStrategy, VcfRecord, VcfStore};
use crossbeam_channel::{bounded, Receiver, Sender};
use rayon::prelude::*;
use std::collections::{HashMap, VecDeque};
use std::path::PathBuf;
use std::rc::Rc;
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

struct SameLenBatchPatch<'a> {
    idx: usize,
    rlen: usize,
    ref_allele: &'a [u8],
    alt: &'a [u8],
    gt_bits: &'a crate::vcf_store::BiallelicPhasedGtBits,
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
    opts: EngineOptions,
}

impl ConsensusEngine {
    pub fn new(ref_index: RefIndex, vcfs: HashMap<String, VcfStore>, opts: EngineOptions) -> Self {
        let vcfs = vcfs.into_iter().map(|(k, v)| (k, Arc::new(v))).collect();
        ConsensusEngine {
            ref_index: Arc::new(ref_index),
            vcfs: Arc::new(vcfs),
            opts,
        }
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
        Ok(ConsensusEngine::new(ref_index, vcfs, opts))
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
        let plan_opts = self.plan_options();
        let (records, plan) = vcf.plan_query(
            &task.chr,
            ori_pos,
            end0,
            self.opts.regions_overlap,
            plan_opts,
        );

        // Build sample mode + options for this task.
        let sample_mode =
            build_sample_mode(vcf, &task.sample, &task.haplotype, self.opts.iupac_codes);
        let mask = match &self.opts.mask {
            Some(p) => match crate::mask::Mask::load(p, self.opts.mask_with) {
                Ok(m) => Some(Rc::new(m)),
                Err(e) => return mk_err(task, format!("mask load failed: {}", e)),
            },
            None => None,
        };
        let opts = ApplyOptions {
            absent_allele: self.opts.absent,
            missing_allele: self.opts.missing,
            mark_del: self.opts.mark_del,
            mark_ins: self.opts.mark_ins,
            mark_snv: self.opts.mark_snv,
            sample_mode,
            mask,
        };

        if self.opts.chain {
            let mut chain = Chain::new(task.chr.clone(), ori_pos, ref_seq.len() as i64);
            let state = apply_region_planned(
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
            let state = apply_region_planned(
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
        let plan_opts = self.plan_options();
        let (records, plan) = vcf.plan_query(
            &group.key.chr,
            ori_pos,
            end0,
            self.opts.regions_overlap,
            plan_opts,
        );
        if let Some(batch) = self
            .try_run_biallelic_phased_batch(tasks, group, vcf, &ref_seq, ori_pos, &records, &plan)
        {
            return batch;
        }
        let shared_mask = match &self.opts.mask {
            Some(p) => match crate::mask::Mask::load(p, self.opts.mask_with) {
                Ok(m) => Some(Rc::new(m)),
                Err(e) => return group_error(format!("mask load failed: {}", e)),
            },
            None => None,
        };

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

            let key = task_exec_key(task);
            if cache_enabled {
                if let Some(cached) = output_cache.get(&key) {
                    out.push((idx, result_from_cached(task, cached)));
                    continue;
                }
            }

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
            if cache_enabled {
                output_cache.insert(key, CachedOutput::from(&result));
            }
            out.push((idx, result));
        }
        out
    }

    fn run_group_task(
        &self,
        task: &ConsensusTask,
        vcf: &VcfStore,
        ref_for_task: Vec<u8>,
        ori_pos: i64,
        records: &[&VcfRecord],
        plan: &RegionPlan,
        shared_mask: Option<Rc<crate::mask::Mask>>,
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
            let state = apply_region_planned(
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
            let state = apply_region_planned(
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

    fn try_run_biallelic_phased_batch(
        &self,
        tasks: &[ConsensusTask],
        group: &TaskGroup,
        vcf: &VcfStore,
        ref_seq: &[u8],
        ori_pos: i64,
        records: &[&VcfRecord],
        plan: &RegionPlan,
    ) -> Option<Vec<(usize, ConsensusResult)>> {
        if records.is_empty()
            || plan.lane != FastPathLane::SameLenOnly
            || self.opts.chain
            || self.opts.mask.is_some()
            || self.opts.absent.is_some()
            || self.opts.missing.is_some()
        {
            return None;
        }

        let mut batch_tasks = Vec::with_capacity(group.indices.len());
        let mut batch_index: HashMap<(usize, usize), usize> = HashMap::new();
        let mut output_map = Vec::with_capacity(group.indices.len());
        for &input_idx in &group.indices {
            let task = &tasks[input_idx];
            let sample = task.sample.as_ref()?;
            let sample_idx = vcf.sample_index(sample)?;
            if sample_idx < 0 {
                return None;
            }
            let hap = task
                .haplotype
                .as_ref()
                .and_then(|h| HaplotypeSpec::parse(h))
                .and_then(|spec| spec.haplotype)?;
            if hap == 0 || hap > 2 {
                return None;
            }
            let batch_key = (sample_idx as usize, hap as usize);
            let batch_idx = match batch_index.get(&batch_key) {
                Some(&idx) => idx,
                None => {
                    let idx = batch_tasks.len();
                    batch_index.insert(batch_key, idx);
                    batch_tasks.push(BatchTask {
                        sample_idx: batch_key.0,
                        hap: batch_key.1,
                    });
                    idx
                }
            };
            output_map.push((input_idx, batch_idx));
        }
        if batch_tasks.is_empty() {
            return None;
        }

        let mut patches = Vec::with_capacity(records.len());
        let mut frz_pos = -1i64;
        for rec in records {
            if rec.alleles.len() == 1 {
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
            if idx + rlen > ref_seq.len() {
                return None;
            }
            if !ref_seq[idx..idx + rlen].eq_ignore_ascii_case(&rec.alleles[0]) {
                return None;
            }
            patches.push(SameLenBatchPatch {
                idx,
                rlen,
                ref_allele: &rec.alleles[0],
                alt: &rec.alleles[1],
                gt_bits,
            });
            frz_pos = rec.ref_end();
        }
        if patches.is_empty() {
            return None;
        }

        let mut buffers: Vec<Option<Vec<u8>>> =
            batch_tasks.iter().map(|_| Some(ref_seq.to_vec())).collect();
        for patch in &patches {
            let to_upper = ref_seq[patch.idx].is_ascii_uppercase();
            for (task_i, task) in batch_tasks.iter().enumerate() {
                if !patch
                    .gt_bits
                    .is_alt_for_hap(task.sample_idx, task.hap)
                    .unwrap_or(false)
                {
                    continue;
                }
                let buf = buffers[task_i].as_mut().expect("batch buffer present");
                let dst = &mut buf[patch.idx..patch.idx + patch.rlen];
                copy_alt_with_case(dst, patch.alt, to_upper);
                if let Some(mark) = self.opts.mark_snv {
                    mark_snv_in_place(patch.ref_allele, dst, mark);
                }
            }
        }

        let mut remaining = vec![0usize; batch_tasks.len()];
        for &(_, batch_idx) in &output_map {
            remaining[batch_idx] += 1;
        }
        let mut out = Vec::with_capacity(output_map.len());
        for (input_idx, batch_idx) in output_map {
            let src_task = &tasks[input_idx];
            remaining[batch_idx] -= 1;
            let seq = if remaining[batch_idx] == 0 {
                buffers[batch_idx]
                    .take()
                    .expect("last use consumes batch buffer")
            } else {
                buffers[batch_idx]
                    .as_ref()
                    .expect("batch buffer available")
                    .clone()
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

    fn plan_options(&self) -> crate::planner::PlanOptions {
        crate::planner::PlanOptions {
            chain: self.opts.chain,
            mark_del: self.opts.mark_del.is_some(),
            mark_ins: self.opts.mark_ins.is_some(),
            mark_snv: self.opts.mark_snv.is_some(),
            mask: self.opts.mask.is_some(),
            absent: self.opts.absent.is_some(),
        }
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
            let _ = HaplotypeSpec::parse(h);
            SampleMode::ApplyAllAlt
        }
    }
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
    let mut seen: HashMap<TaskExecKey, ()> = HashMap::with_capacity(indices.len());
    for &idx in indices {
        let key = task_exec_key(&tasks[idx]);
        if seen.insert(key, ()).is_some() {
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
        alt.iter().any(u8::is_ascii_lowercase)
    } else {
        alt.iter().any(u8::is_ascii_uppercase)
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
        let (records, plan) = vcf.plan_query(
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
        let (records, plan) = vcf.plan_query(
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
}
