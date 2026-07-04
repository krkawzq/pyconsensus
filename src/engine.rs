//! engine — ConsensusEngine: holds preprocessed material (ref + vcfs), runs
//! multi-threaded consensus production.
//!
//! (docs/design.md §5.4 / §5.6) Each (region, sample, haplotype) task is
//! independent → natural parallelism. `consensus_many` and `consensus_iter`
//! group identical regions so ref fetch, VCF query, and region planning are
//! amortized across sample/haplotype tasks.
//!
//! This module is PyO3-free; `py.rs` wraps it under the `python` feature.

use crate::apply::{
    apply_region_planned_set, apply_region_planned_set_profile, apply_region_planned_slice_profile,
    force_fallback_state_machine, ApplyOptions, TO_LOWER, TO_UPPER,
};
use crate::chain::Chain;
use crate::compiled::{
    allele_case_flags, RecordFlags, ALLELE_HAS_ASCII_LOWER, ALLELE_HAS_ASCII_UPPER,
};
use crate::haplotype::{HaplotypeSpec, SampleMode};
use crate::planner::{plan_region_set, PlanOptions, RegionPlan};
use crate::ref_index::RefIndex;
use crate::stats::{FallbackReason, FastPathLane, RuntimeStats};
use crate::vcf_store::{LoadStrategy, RecordSet, VcfStore};
use crossbeam_channel::{bounded, Receiver, Sender};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet, VecDeque};
use std::path::PathBuf;
use std::sync::{Arc, OnceLock};

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
    missing_out: u8,
    gt_bits: &'a crate::vcf_store::BiallelicPhasedGtBits,
}

struct MnpBatchPatch<'a> {
    idx: usize,
    rlen: usize,
    ref_allele: &'a [u8],
    ref_case_flags: u8,
    alt: &'a [u8],
    alt_case_flags: u8,
    missing_out: u8,
    gt_bits: &'a crate::vcf_store::BiallelicPhasedGtBits,
    to_upper: bool,
}

enum SameLenBatchPatch<'a> {
    RefOnly(RefOnlyBatchPatch<'a>),
    Snp1(Snp1BatchPatch<'a>),
    Mnp(MnpBatchPatch<'a>),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BatchExecutionFlavor {
    Plain,
    Missing { missing: u8 },
    Absent { absent: u8 },
    AbsentMissing { absent: u8, missing: u8 },
}

#[inline]
fn biallelic_batch_fastpath_disabled() -> bool {
    static DISABLED: OnceLock<bool> = OnceLock::new();
    force_fallback_state_machine()
        || *DISABLED.get_or_init(|| {
            std::env::var_os("PYCONSENSUS_DISABLE_BIALLELIC_BATCH_FASTPATH").is_some()
        })
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

/// Aggregate over produced consensus results for blackhole/throughput runs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConsensusRunStats {
    pub n: usize,
    pub total_len: u64,
    pub min_len: usize,
    pub max_len: usize,
}

impl Default for ConsensusRunStats {
    fn default() -> Self {
        ConsensusRunStats {
            n: 0,
            total_len: 0,
            min_len: usize::MAX,
            max_len: 0,
        }
    }
}

impl ConsensusRunStats {
    fn observe(&mut self, result: ConsensusResult) -> Result<(), String> {
        if let Some(err) = result.error {
            return Err(format!("{}: {}", result.gene_id, err));
        }
        let len = result.seq.len();
        self.n += 1;
        self.total_len += len as u64;
        self.min_len = self.min_len.min(len);
        self.max_len = self.max_len.max(len);
        Ok(())
    }

    fn merge(&mut self, other: ConsensusRunStats) {
        if other.n == 0 {
            return;
        }
        self.n += other.n;
        self.total_len += other.total_len;
        self.min_len = self.min_len.min(other.min_len);
        self.max_len = self.max_len.max(other.max_len);
    }

    fn finish(mut self) -> Self {
        if self.n == 0 {
            self.min_len = 0;
        }
        self
    }

    pub fn as_tuple(&self) -> (usize, u64, usize, usize) {
        (self.n, self.total_len, self.min_len, self.max_len)
    }
}

/// Combined throughput and dispatch counters for a consensus run.
#[derive(Clone, Debug, Default)]
pub struct ConsensusRunProfile {
    pub run: ConsensusRunStats,
    pub runtime: RuntimeStats,
    pub elapsed_secs: f64,
}

impl ConsensusRunProfile {
    fn observe_result(&mut self, result: ConsensusResult) -> Result<(), String> {
        self.runtime
            .observe_alloc_bytes(u64::try_from(result.seq.len()).unwrap_or(u64::MAX));
        self.run.observe(result)
    }

    fn merge(&mut self, other: ConsensusRunProfile) {
        self.run.merge(other.run);
        self.runtime.merge(other.runtime);
    }

    fn finish(mut self) -> Self {
        self.run = self.run.finish();
        self
    }

    pub fn summary_lines(&self) -> Vec<String> {
        let mut lines = vec![
            format!("run.n={}", self.run.n),
            format!("run.total_len={}", self.run.total_len),
            format!("run.min_len={}", self.run.min_len),
            format!("run.max_len={}", self.run.max_len),
            format!("run.elapsed_sec={:.6}", self.elapsed_secs),
            format!(
                "run.seq_per_sec={:.3}",
                rate(self.run.n as f64, self.elapsed_secs)
            ),
            format!(
                "run.bases_per_sec={:.3}",
                rate(self.run.total_len as f64, self.elapsed_secs)
            ),
            format!(
                "runtime.records_per_sec={:.3}",
                rate(self.runtime.records_seen as f64, self.elapsed_secs)
            ),
            format!(
                "runtime.tasks_per_sec={:.3}",
                rate(self.runtime.tasks_total as f64, self.elapsed_secs)
            ),
        ];
        lines.extend(
            self.runtime
                .summary_lines()
                .into_iter()
                .map(|line| format!("runtime.{}", line)),
        );
        lines
    }
}

struct ProfiledTaskResult {
    result: ConsensusResult,
    lane: FastPathLane,
    napplied: u64,
    fallback_reason: Option<FallbackReason>,
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
    /// 0 = unlimited; otherwise split each region/VCF group into chunks of at
    /// most this many tasks.
    pub max_tasks_per_group: usize,
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

        let groups = group_tasks(&tasks, self.opts.max_tasks_per_group);
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

    /// Run all tasks in parallel and consume results inside Rust.
    ///
    /// This is intended for throughput tests and sinks where callers do not
    /// need every sequence materialized as Python objects. It keeps only a
    /// small aggregate per worker group instead of the full ordered result set.
    pub fn consensus_many_stats(
        &self,
        tasks: Vec<ConsensusTask>,
        threads: usize,
    ) -> Result<ConsensusRunStats, String> {
        if tasks.is_empty() {
            return Ok(ConsensusRunStats::default().finish());
        }
        let nthr = threads.max(1);
        let pool = thread_pool(nthr);

        let groups = group_tasks(&tasks, self.opts.max_tasks_per_group);
        let tasks = Arc::new(tasks);
        let engine = self.clone_shallow();
        let stats = pool.install(move || {
            groups
                .par_iter()
                .map(|group| {
                    let mut stats = ConsensusRunStats::default();
                    for (_, result) in engine.run_group(&tasks, group) {
                        stats.observe(result)?;
                    }
                    Ok(stats)
                })
                .reduce(|| Ok(ConsensusRunStats::default()), merge_run_stats)
        })?;
        Ok(stats.finish())
    }

    /// Run tasks and return dispatch counters without materializing results to
    /// the caller. This is the observability path for benchmark harnesses:
    /// region grouping and biallelic batch execution remain enabled, while the
    /// returned profile records actual lanes used by each group/task.
    pub fn consensus_many_profile(
        &self,
        tasks: Vec<ConsensusTask>,
        threads: usize,
    ) -> Result<ConsensusRunProfile, String> {
        if tasks.is_empty() {
            return Ok(ConsensusRunProfile::default().finish());
        }
        let started = std::time::Instant::now();
        let nthr = threads.max(1);
        let pool = thread_pool(nthr);

        let groups = group_tasks(&tasks, self.opts.max_tasks_per_group);
        let tasks = Arc::new(tasks);
        let engine = self.clone_shallow();
        let mut profile = pool.install(move || {
            groups
                .par_iter()
                .map(|group| engine.run_group_profile(&tasks, group))
                .reduce(|| Ok(ConsensusRunProfile::default()), merge_run_profiles)
        })?;
        profile.elapsed_secs = started.elapsed().as_secs_f64();
        Ok(profile.finish())
    }

    /// VCF-load compile counters, one key-value line per counter. The key is
    /// prefixed by `vcf.<vcf_key>.` so callers can log the vector directly.
    pub fn compile_stats_lines(&self) -> Vec<String> {
        let mut keys: Vec<_> = self.vcfs.keys().cloned().collect();
        keys.sort();
        let mut lines = Vec::new();
        for key in keys {
            if let Some(vcf) = self.vcfs.get(&key) {
                lines.extend(
                    vcf.compile_stats()
                        .summary_lines()
                        .into_iter()
                        .map(|line| format!("vcf.{}.{}", key, line)),
                );
            }
        }
        lines
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
        let n = group.indices.len();
        let mut out = Vec::with_capacity(n);
        let use_borrowed_ref = can_use_borrowed_ref(&plan, self.opts.chain);
        let mut owned_ref_seq = Some(ref_seq);
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

            let result = if use_borrowed_ref {
                self.run_group_task_borrowed(
                    task,
                    vcf,
                    owned_ref_seq
                        .as_deref()
                        .expect("borrowed ref path keeps shared ref"),
                    ori_pos,
                    &records,
                    &plan,
                    shared_mask.clone(),
                )
            } else {
                let ref_for_task = if cache_enabled {
                    owned_ref_seq
                        .as_ref()
                        .expect("shared ref available")
                        .clone()
                } else if j + 1 == n {
                    owned_ref_seq.take().expect("last task consumes shared ref")
                } else {
                    owned_ref_seq
                        .as_ref()
                        .expect("shared ref available")
                        .clone()
                };
                self.run_group_task(
                    task,
                    vcf,
                    ref_for_task,
                    ori_pos,
                    &records,
                    &plan,
                    shared_mask.clone(),
                )
            };
            if let Some(key) = cache_key {
                output_cache.insert(key, CachedOutput::from(&result));
            }
            out.push((idx, result));
        }
        out
    }

    fn run_group_profile(
        &self,
        tasks: &[ConsensusTask],
        group: &TaskGroup,
    ) -> Result<ConsensusRunProfile, String> {
        let mut profile = ConsensusRunProfile::default();
        profile.runtime.observe_region();
        profile
            .runtime
            .observe_tasks(u64::try_from(group.indices.len()).unwrap_or(u64::MAX));

        let group_error = |err: String| -> String {
            let first = &tasks[group.indices[0]];
            format!(
                "{}:{}-{}:{}: {}",
                first.chr, first.start, first.end, first.vcf_key, err
            )
        };

        let vcf = self
            .vcfs
            .get(&group.key.vcf_key)
            .ok_or_else(|| group_error(format!("unknown vcf_key: {}", group.key.vcf_key)))?;

        let ref_seq = self
            .ref_index
            .fetch_1based(&group.key.chr, group.key.start, group.key.end)
            .map_err(|e| group_error(format!("ref fetch failed: {}", e)))?;
        let ori_pos = group.key.start - 1;
        let end0 = group.key.end - 1;
        let records = vcf.query_set(&group.key.chr, ori_pos, end0, self.opts.regions_overlap);
        let plan_opts = self.plan_options_for_records(&group.key.chr, &records);
        let plan = plan_region_set(&records, plan_opts);

        if let Some(batch) = self
            .try_run_biallelic_phased_batch(tasks, group, vcf, &ref_seq, ori_pos, &records, &plan)
        {
            profile
                .runtime
                .observe_records(u64::try_from(records.len()).unwrap_or(u64::MAX));
            profile
                .runtime
                .observe_lane(FastPathLane::BiallelicPhasedBatch);
            profile.runtime.observe_same_len_fastpath_records(
                u64::try_from(records.len()).unwrap_or(u64::MAX),
            );
            for (_, result) in batch {
                profile.observe_result(result)?;
            }
            return Ok(profile);
        }

        let shared_mask = self.mask.clone();
        if records.is_empty() {
            profile.runtime.observe_lane(FastPathLane::EmptyRegion);
            for (_, result) in
                self.run_empty_group(tasks, group, &ref_seq, ori_pos, shared_mask.as_deref())
            {
                profile.observe_result(result)?;
            }
            return Ok(profile);
        }

        let cache_enabled = has_duplicate_exec_keys(tasks, &group.indices);
        let mut output_cache: HashMap<TaskExecKey, CachedOutput> = HashMap::new();
        let n = group.indices.len();
        let use_borrowed_ref = can_use_borrowed_ref(&plan, self.opts.chain);
        let mut owned_ref_seq = Some(ref_seq);
        for (j, &idx) in group.indices.iter().enumerate() {
            let task = &tasks[idx];
            let cache_key = if cache_enabled {
                let key = task_exec_key(task);
                if let Some(cached) = output_cache.get(&key) {
                    profile.observe_result(result_from_cached(task, cached))?;
                    continue;
                }
                Some(key)
            } else {
                None
            };

            profile
                .runtime
                .observe_records(u64::try_from(records.len()).unwrap_or(u64::MAX));
            let observed = if use_borrowed_ref {
                self.run_group_task_borrowed_profile(
                    task,
                    vcf,
                    owned_ref_seq
                        .as_deref()
                        .expect("borrowed ref path keeps shared ref"),
                    ori_pos,
                    &records,
                    &plan,
                    shared_mask.clone(),
                )
            } else {
                let ref_for_task = if cache_enabled {
                    owned_ref_seq
                        .as_ref()
                        .expect("shared ref available")
                        .clone()
                } else if j + 1 == n {
                    owned_ref_seq.take().expect("last task consumes shared ref")
                } else {
                    owned_ref_seq
                        .as_ref()
                        .expect("shared ref available")
                        .clone()
                };
                self.run_group_task_profile(
                    task,
                    vcf,
                    ref_for_task,
                    ori_pos,
                    &records,
                    &plan,
                    shared_mask.clone(),
                )
            };
            observe_profiled_task(
                &mut profile,
                observed.lane,
                observed.napplied,
                observed.fallback_reason,
                records.len(),
                &plan,
            );
            if let Some(key) = cache_key {
                output_cache.insert(key, CachedOutput::from(&observed.result));
            }
            profile.observe_result(observed.result)?;
        }
        Ok(profile)
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
    fn run_group_task_borrowed(
        &self,
        task: &ConsensusTask,
        vcf: &VcfStore,
        ref_seq: &[u8],
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

        let (state, chain) = if self.opts.chain {
            let mut chain = Chain::new(task.chr.clone(), ori_pos, ref_seq.len() as i64);
            let (state, _lane) = apply_region_planned_slice_profile(
                &task.chr,
                ref_seq,
                ori_pos,
                records,
                &opts,
                Some(&mut chain),
                Some(plan),
            );
            (state, Some(chain.render()))
        } else {
            let (state, _lane) = apply_region_planned_slice_profile(
                &task.chr,
                ref_seq,
                ori_pos,
                records,
                &opts,
                None,
                Some(plan),
            );
            (state, None)
        };
        ConsensusResult {
            gene_id: task.gene_id.clone(),
            sample: task.sample.clone(),
            haplotype: task.haplotype.clone(),
            seq: state.buf,
            chain,
            error: None,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn run_group_task_profile(
        &self,
        task: &ConsensusTask,
        vcf: &VcfStore,
        ref_for_task: Vec<u8>,
        ori_pos: i64,
        records: &RecordSet<'_>,
        plan: &RegionPlan,
        shared_mask: Option<Arc<crate::mask::Mask>>,
    ) -> ProfiledTaskResult {
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

        let (state, lane, chain) = if self.opts.chain {
            let mut chain = Chain::new(task.chr.clone(), ori_pos, ref_for_task.len() as i64);
            let (state, lane) = apply_region_planned_set_profile(
                &task.chr,
                ref_for_task,
                ori_pos,
                records,
                &opts,
                Some(&mut chain),
                Some(plan),
            );
            (state, lane, Some(chain.render()))
        } else {
            let (state, lane) = apply_region_planned_set_profile(
                &task.chr,
                ref_for_task,
                ori_pos,
                records,
                &opts,
                None,
                Some(plan),
            );
            (state, lane, None)
        };
        let napplied = state.napplied;
        let fallback_reason = state.fallback_reason;
        ProfiledTaskResult {
            result: ConsensusResult {
                gene_id: task.gene_id.clone(),
                sample: task.sample.clone(),
                haplotype: task.haplotype.clone(),
                seq: state.buf,
                chain,
                error: None,
            },
            lane,
            napplied,
            fallback_reason,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn run_group_task_borrowed_profile(
        &self,
        task: &ConsensusTask,
        vcf: &VcfStore,
        ref_seq: &[u8],
        ori_pos: i64,
        records: &RecordSet<'_>,
        plan: &RegionPlan,
        shared_mask: Option<Arc<crate::mask::Mask>>,
    ) -> ProfiledTaskResult {
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

        let (state, lane, chain) = if self.opts.chain {
            let mut chain = Chain::new(task.chr.clone(), ori_pos, ref_seq.len() as i64);
            let (state, lane) = apply_region_planned_slice_profile(
                &task.chr,
                ref_seq,
                ori_pos,
                records,
                &opts,
                Some(&mut chain),
                Some(plan),
            );
            (state, lane, Some(chain.render()))
        } else {
            let (state, lane) = apply_region_planned_slice_profile(
                &task.chr,
                ref_seq,
                ori_pos,
                records,
                &opts,
                None,
                Some(plan),
            );
            (state, lane, None)
        };
        let napplied = state.napplied;
        let fallback_reason = state.fallback_reason;
        ProfiledTaskResult {
            result: ConsensusResult {
                gene_id: task.gene_id.clone(),
                sample: task.sample.clone(),
                haplotype: task.haplotype.clone(),
                seq: state.buf,
                chain,
                error: None,
            },
            lane,
            napplied,
            fallback_reason,
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
        if biallelic_batch_fastpath_disabled()
            || records.is_empty()
            || plan.lane != FastPathLane::SameLenOnly
            || self.opts.chain
        {
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
        let mut n_batch_outputs = 0usize;
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
                let idx = n_batch_outputs;
                *slot = idx;
                n_batch_outputs += 1;
                idx
            } else {
                *slot
            };
            output_map.push((input_idx, batch_idx));
        }
        if n_batch_outputs == 0 {
            return None;
        }

        let mut patches = Vec::with_capacity(records.len());
        let mut frz_pos = -1i64;
        for rec_idx in records.iter_indices()? {
            let (pos, rlen_i32, ref_end) = vcf.compiled_span(rec_idx)?;
            let n_alleles = vcf.compiled_n_alleles(rec_idx)?;
            if n_alleles == 1 {
                if self.opts.absent.is_none() {
                    continue;
                }
                if pos <= frz_pos || rlen_i32 <= 0 {
                    return None;
                }
                let rlen = rlen_i32 as usize;
                let ref_allele = vcf.compiled_allele(rec_idx, 0)?;
                if ref_allele.len() != rlen {
                    return None;
                }
                let idx = pos - ori_pos;
                if idx < 0 {
                    return None;
                }
                let idx = idx as usize;
                if idx + rlen > base_ref.len() {
                    return None;
                }
                if !base_ref[idx..idx + rlen].eq_ignore_ascii_case(ref_allele) {
                    return None;
                }
                let ref_case_flags = vcf
                    .compiled_allele_op(rec_idx, 0)
                    .map(|op| op.case_flags)
                    .unwrap_or_else(|| allele_case_flags(ref_allele));
                let to_upper = base_ref[idx].is_ascii_uppercase();
                patches.push(SameLenBatchPatch::RefOnly(RefOnlyBatchPatch {
                    idx,
                    rlen,
                    ref_allele,
                    ref_case_flags,
                    ref_out: snp1_alt_with_case_and_mark(
                        ref_allele[0],
                        ref_allele[0],
                        to_upper,
                        ref_case_flags,
                        None,
                    ),
                    to_upper,
                }));
                frz_pos = ref_end;
                continue;
            }
            let flags = vcf.compiled_flags(rec_idx)?;
            if pos <= frz_pos
                || rlen_i32 <= 0
                || !flags.contains(RecordFlags::BIALLELIC)
                || !flags.contains(RecordFlags::ALL_ALT_SAME_LEN)
            {
                return None;
            }
            let gt_bits = vcf.compiled_gt_bits(rec_idx)?;
            if active_samples_need_gt_fallback(
                gt_bits,
                &active_words_by_hap,
                &active_word_indices_by_hap,
            ) {
                return None;
            }
            let rlen = rlen_i32 as usize;
            let ref_allele = vcf.compiled_allele(rec_idx, 0)?;
            let alt_allele = vcf.compiled_allele(rec_idx, 1)?;
            if n_alleles != 2 || ref_allele.len() != rlen || alt_allele.len() != rlen {
                return None;
            }
            let idx = pos - ori_pos;
            if idx < 0 {
                return None;
            }
            let idx = idx as usize;
            if idx + rlen > base_ref.len() {
                return None;
            }
            if !base_ref[idx..idx + rlen].eq_ignore_ascii_case(ref_allele) {
                return None;
            }
            let ref_case_flags = vcf
                .compiled_allele_op(rec_idx, 0)
                .map(|op| op.case_flags)
                .unwrap_or_else(|| allele_case_flags(ref_allele));
            let alt_case_flags = vcf
                .compiled_allele_op(rec_idx, 1)
                .map(|op| op.case_flags)
                .unwrap_or_else(|| allele_case_flags(alt_allele));
            let to_upper = base_ref[idx].is_ascii_uppercase();
            if rlen == 1 {
                let ref_base = ref_allele[0];
                let alt_base = alt_allele[0];
                let missing_out = self
                    .opts
                    .missing
                    .map(|missing| {
                        snp1_alt_with_case_and_mark(
                            ref_base,
                            missing,
                            to_upper,
                            byte_case_flags(missing),
                            self.opts.mark_snv,
                        )
                    })
                    .unwrap_or(0);
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
                    missing_out,
                    gt_bits,
                }));
            } else {
                let ref_base = ref_allele[0];
                let missing_out = self
                    .opts
                    .missing
                    .map(|missing| {
                        snp1_alt_with_case_and_mark(
                            ref_base,
                            missing,
                            to_upper,
                            byte_case_flags(missing),
                            self.opts.mark_snv,
                        )
                    })
                    .unwrap_or(0);
                patches.push(SameLenBatchPatch::Mnp(MnpBatchPatch {
                    idx,
                    rlen,
                    ref_allele,
                    ref_case_flags,
                    alt: alt_allele,
                    alt_case_flags,
                    missing_out,
                    gt_bits,
                    to_upper,
                }));
            }
            frz_pos = ref_end;
        }
        if patches.is_empty() {
            return None;
        }

        match self.biallelic_batch_flavor() {
            BatchExecutionFlavor::Plain => self.try_run_biallelic_phased_alt_only_batch(
                tasks,
                base_ref,
                &output_map,
                &task_by_hap_sample,
                &active_words_by_hap,
                &active_word_indices_by_hap,
                &patches,
            ),
            BatchExecutionFlavor::Missing { missing } => self
                .try_run_biallelic_phased_missing_batch(
                    tasks,
                    base_ref,
                    &output_map,
                    &task_by_hap_sample,
                    &active_words_by_hap,
                    &active_word_indices_by_hap,
                    &patches,
                    missing,
                ),
            BatchExecutionFlavor::Absent { absent } => self.try_run_biallelic_phased_absent_batch(
                tasks,
                base_ref,
                &output_map,
                &task_by_hap_sample,
                &active_words_by_hap,
                &active_word_indices_by_hap,
                &patches,
                absent,
            ),
            BatchExecutionFlavor::AbsentMissing { absent, missing } => self
                .try_run_biallelic_phased_absent_missing_batch(
                    tasks,
                    base_ref,
                    &output_map,
                    &task_by_hap_sample,
                    &active_words_by_hap,
                    &active_word_indices_by_hap,
                    &patches,
                    absent,
                    missing,
                ),
        }
    }

    fn biallelic_batch_flavor(&self) -> BatchExecutionFlavor {
        match (self.opts.absent, self.opts.missing) {
            (None, None) => BatchExecutionFlavor::Plain,
            (None, Some(missing)) => BatchExecutionFlavor::Missing { missing },
            (Some(absent), None) => BatchExecutionFlavor::Absent { absent },
            (Some(absent), Some(missing)) => {
                BatchExecutionFlavor::AbsentMissing { absent, missing }
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
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
        if output_map.is_empty() {
            return None;
        }
        let n_buffers = batch_buffer_count(output_map)?;
        let mut buffers: Vec<Vec<u8>> = (0..n_buffers).map(|_| base_ref.to_vec()).collect();
        for patch in patches {
            match patch {
                SameLenBatchPatch::RefOnly(_) => continue,
                SameLenBatchPatch::Snp1(patch) => {
                    for hap_idx in 0..2 {
                        let task_by_sample = &task_by_hap_sample[hap_idx];
                        let active_words = &active_words_by_hap[hap_idx];
                        let active_word_indices = &active_word_indices_by_hap[hap_idx];
                        let words = patch.gt_bits.alt_words_for_hap_index(hap_idx);
                        for &word_idx in active_word_indices {
                            let mut bits = words[word_idx] & active_words[word_idx];
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
                    for hap_idx in 0..2 {
                        let task_by_sample = &task_by_hap_sample[hap_idx];
                        let active_words = &active_words_by_hap[hap_idx];
                        let active_word_indices = &active_word_indices_by_hap[hap_idx];
                        let words = patch.gt_bits.alt_words_for_hap_index(hap_idx);
                        for &word_idx in active_word_indices {
                            let mut bits = words[word_idx] & active_words[word_idx];
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

        Some(finish_batch_outputs(tasks, output_map, buffers))
    }

    #[allow(clippy::too_many_arguments)]
    fn try_run_biallelic_phased_missing_batch(
        &self,
        tasks: &[ConsensusTask],
        base_ref: &[u8],
        output_map: &[(usize, usize)],
        task_by_hap_sample: &[Vec<usize>; 2],
        active_words_by_hap: &[Vec<u64>; 2],
        active_word_indices_by_hap: &[Vec<usize>; 2],
        patches: &[SameLenBatchPatch<'_>],
        _missing: u8,
    ) -> Option<Vec<(usize, ConsensusResult)>> {
        if output_map.is_empty() {
            return None;
        }
        let n_buffers = batch_buffer_count(output_map)?;
        let mut buffers: Vec<Vec<u8>> = (0..n_buffers).map(|_| base_ref.to_vec()).collect();
        for patch in patches {
            match patch {
                SameLenBatchPatch::RefOnly(_) => continue,
                SameLenBatchPatch::Snp1(patch) => {
                    for hap_idx in 0..2 {
                        let task_by_sample = &task_by_hap_sample[hap_idx];
                        let active_words = &active_words_by_hap[hap_idx];
                        let active_word_indices = &active_word_indices_by_hap[hap_idx];
                        let alt_words = patch.gt_bits.alt_words_for_hap_index(hap_idx);
                        let missing_words = patch.gt_bits.missing_words_for_hap_index(hap_idx);
                        for &word_idx in active_word_indices {
                            let active = active_words[word_idx];
                            let mut alt_bits = alt_words[word_idx] & active;
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

                            let mut missing_bits = missing_words[word_idx] & active;
                            while missing_bits != 0 {
                                let bit_idx = missing_bits.trailing_zeros() as usize;
                                let sample_idx = word_idx * 64 + bit_idx;
                                let task_idx = task_by_sample[sample_idx];
                                debug_assert_ne!(task_idx, usize::MAX);
                                if task_idx != usize::MAX {
                                    buffers[task_idx][patch.idx] = patch.missing_out;
                                }
                                missing_bits &= missing_bits - 1;
                            }
                        }
                    }
                }
                SameLenBatchPatch::Mnp(patch) => {
                    for hap_idx in 0..2 {
                        let task_by_sample = &task_by_hap_sample[hap_idx];
                        let active_words = &active_words_by_hap[hap_idx];
                        let active_word_indices = &active_word_indices_by_hap[hap_idx];
                        let alt_words = patch.gt_bits.alt_words_for_hap_index(hap_idx);
                        let missing_words = patch.gt_bits.missing_words_for_hap_index(hap_idx);
                        for &word_idx in active_word_indices {
                            let active = active_words[word_idx];
                            let mut alt_bits = alt_words[word_idx] & active;
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

                            let mut missing_bits = missing_words[word_idx] & active;
                            while missing_bits != 0 {
                                let bit_idx = missing_bits.trailing_zeros() as usize;
                                let sample_idx = word_idx * 64 + bit_idx;
                                let task_idx = task_by_sample[sample_idx];
                                debug_assert_ne!(task_idx, usize::MAX);
                                if task_idx != usize::MAX {
                                    let buf = &mut buffers[task_idx];
                                    buf[patch.idx] = patch.missing_out;
                                }
                                missing_bits &= missing_bits - 1;
                            }
                        }
                    }
                }
            }
        }

        Some(finish_batch_outputs(tasks, output_map, buffers))
    }

    #[allow(clippy::too_many_arguments)]
    fn try_run_biallelic_phased_absent_batch(
        &self,
        tasks: &[ConsensusTask],
        base_ref: &[u8],
        output_map: &[(usize, usize)],
        task_by_hap_sample: &[Vec<usize>; 2],
        active_words_by_hap: &[Vec<u64>; 2],
        active_word_indices_by_hap: &[Vec<usize>; 2],
        patches: &[SameLenBatchPatch<'_>],
        absent: u8,
    ) -> Option<Vec<(usize, ConsensusResult)>> {
        if output_map.is_empty() {
            return None;
        }
        let n_buffers = batch_buffer_count(output_map)?;
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
                    for hap_idx in 0..2 {
                        let task_by_sample = &task_by_hap_sample[hap_idx];
                        let active_words = &active_words_by_hap[hap_idx];
                        let active_word_indices = &active_word_indices_by_hap[hap_idx];
                        let alt_words = patch.gt_bits.alt_words_for_hap_index(hap_idx);
                        let missing_words = patch.gt_bits.missing_words_for_hap_index(hap_idx);
                        for &word_idx in active_word_indices {
                            let active = active_words[word_idx];
                            let alt_bits = alt_words[word_idx] & active;
                            let missing_bits = missing_words[word_idx] & active;

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
                        }
                    }
                }
                SameLenBatchPatch::Mnp(patch) => {
                    for hap_idx in 0..2 {
                        let task_by_sample = &task_by_hap_sample[hap_idx];
                        let active_words = &active_words_by_hap[hap_idx];
                        let active_word_indices = &active_word_indices_by_hap[hap_idx];
                        let alt_words = patch.gt_bits.alt_words_for_hap_index(hap_idx);
                        let missing_words = patch.gt_bits.missing_words_for_hap_index(hap_idx);
                        for &word_idx in active_word_indices {
                            let active = active_words[word_idx];
                            let alt_bits = alt_words[word_idx] & active;
                            let missing_bits = missing_words[word_idx] & active;

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
                        }
                    }
                }
            }
        }

        Some(finish_batch_outputs(tasks, output_map, buffers))
    }

    #[allow(clippy::too_many_arguments)]
    fn try_run_biallelic_phased_absent_missing_batch(
        &self,
        tasks: &[ConsensusTask],
        base_ref: &[u8],
        output_map: &[(usize, usize)],
        task_by_hap_sample: &[Vec<usize>; 2],
        active_words_by_hap: &[Vec<u64>; 2],
        active_word_indices_by_hap: &[Vec<usize>; 2],
        patches: &[SameLenBatchPatch<'_>],
        absent: u8,
        _missing: u8,
    ) -> Option<Vec<(usize, ConsensusResult)>> {
        if output_map.is_empty() {
            return None;
        }
        let n_buffers = batch_buffer_count(output_map)?;
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
                    for hap_idx in 0..2 {
                        let task_by_sample = &task_by_hap_sample[hap_idx];
                        let active_words = &active_words_by_hap[hap_idx];
                        let active_word_indices = &active_word_indices_by_hap[hap_idx];
                        let alt_words = patch.gt_bits.alt_words_for_hap_index(hap_idx);
                        let missing_words = patch.gt_bits.missing_words_for_hap_index(hap_idx);
                        for &word_idx in active_word_indices {
                            let active = active_words[word_idx];
                            let alt_bits = alt_words[word_idx] & active;
                            let missing_bits = missing_words[word_idx] & active;

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

                            let mut missing_bits = missing_bits;
                            while missing_bits != 0 {
                                let bit_idx = missing_bits.trailing_zeros() as usize;
                                let sample_idx = word_idx * 64 + bit_idx;
                                let task_idx = task_by_sample[sample_idx];
                                debug_assert_ne!(task_idx, usize::MAX);
                                if task_idx != usize::MAX {
                                    buffers[task_idx][patch.idx] = patch.missing_out;
                                }
                                missing_bits &= missing_bits - 1;
                            }
                        }
                    }
                }
                SameLenBatchPatch::Mnp(patch) => {
                    for hap_idx in 0..2 {
                        let task_by_sample = &task_by_hap_sample[hap_idx];
                        let active_words = &active_words_by_hap[hap_idx];
                        let active_word_indices = &active_word_indices_by_hap[hap_idx];
                        let alt_words = patch.gt_bits.alt_words_for_hap_index(hap_idx);
                        let missing_words = patch.gt_bits.missing_words_for_hap_index(hap_idx);
                        for &word_idx in active_word_indices {
                            let active = active_words[word_idx];
                            let alt_bits = alt_words[word_idx] & active;
                            let missing_bits = missing_words[word_idx] & active;

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

                            let mut missing_bits = missing_bits;
                            while missing_bits != 0 {
                                let bit_idx = missing_bits.trailing_zeros() as usize;
                                let sample_idx = word_idx * 64 + bit_idx;
                                let task_idx = task_by_sample[sample_idx];
                                debug_assert_ne!(task_idx, usize::MAX);
                                if task_idx != usize::MAX {
                                    let buf = &mut buffers[task_idx];
                                    buf[patch.idx] = patch.missing_out;
                                }
                                missing_bits &= missing_bits - 1;
                            }
                        }
                    }
                }
            }
        }

        Some(finish_batch_outputs(tasks, output_map, buffers))
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
                    .iter_spans()
                    .any(|span| mask.overlaps(chr, span.pos, span.ref_end))
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
        let groups = Arc::new(group_tasks(&tasks, self.opts.max_tasks_per_group));
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

    /// Drive the lazy iterator to completion while consuming results in Rust.
    pub fn consensus_iter_stats(
        &self,
        tasks: Vec<ConsensusTask>,
        prefetch_steps: usize,
        warmup: bool,
        ordered: bool,
        threads: usize,
    ) -> Result<ConsensusRunStats, String> {
        let mut iter = self.consensus_iter(tasks, prefetch_steps, warmup, ordered, threads);
        let mut stats = ConsensusRunStats::default();
        while let Some((_, result)) = iter.next_blocking() {
            stats.observe(result)?;
        }
        Ok(stats.finish())
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

fn can_use_borrowed_ref(plan: &RegionPlan, chain: bool) -> bool {
    !chain
        && matches!(
            plan.lane,
            FastPathLane::SameLenOnly
                | FastPathLane::SameLenIupac
                | FastPathLane::NormalizedEditScript
                | FastPathLane::MixedSimpleEdits
        )
}

fn merge_run_stats(
    left: Result<ConsensusRunStats, String>,
    right: Result<ConsensusRunStats, String>,
) -> Result<ConsensusRunStats, String> {
    match (left, right) {
        (Ok(mut left), Ok(right)) => {
            left.merge(right);
            Ok(left)
        }
        (Err(err), _) | (_, Err(err)) => Err(err),
    }
}

fn rate(n: f64, elapsed_secs: f64) -> f64 {
    if elapsed_secs > 0.0 {
        n / elapsed_secs
    } else {
        0.0
    }
}

fn merge_run_profiles(
    left: Result<ConsensusRunProfile, String>,
    right: Result<ConsensusRunProfile, String>,
) -> Result<ConsensusRunProfile, String> {
    match (left, right) {
        (Ok(mut left), Ok(right)) => {
            left.merge(right);
            Ok(left)
        }
        (Err(err), _) | (_, Err(err)) => Err(err),
    }
}

fn observe_profiled_task(
    profile: &mut ConsensusRunProfile,
    lane: FastPathLane,
    napplied: u64,
    fallback_reason: Option<FallbackReason>,
    records_len: usize,
    plan: &RegionPlan,
) {
    profile.runtime.observe_lane(lane);
    match lane {
        FastPathLane::SameLenOnly | FastPathLane::SameLenIupac => {
            profile.runtime.observe_same_len_fastpath_records(napplied);
        }
        FastPathLane::BiallelicPhasedBatch => {
            profile
                .runtime
                .observe_same_len_fastpath_records(u64::try_from(records_len).unwrap_or(u64::MAX));
        }
        FastPathLane::NormalizedEditScript | FastPathLane::MixedSimpleEdits => {
            profile
                .runtime
                .observe_edit_script_fastpath_records(napplied);
        }
        FastPathLane::FallbackStateMachine => {
            profile
                .runtime
                .observe_fallback_records(u64::try_from(records_len).unwrap_or(u64::MAX));
            if let Some(reason) = fallback_reason {
                profile.runtime.observe_fallback_reason(reason);
            } else if plan.fallback_reasons.is_empty() {
                profile
                    .runtime
                    .observe_fallback_reason(FallbackReason::UnsupportedMode);
            } else {
                for &reason in &plan.fallback_reasons {
                    profile.runtime.observe_fallback_reason(reason);
                }
            }
        }
        FastPathLane::EmptyRegion => {}
    }
}

struct TaskGroupBuilder {
    key: TaskGroupKey,
    chunks: Vec<Vec<usize>>,
}

fn group_tasks(tasks: &[ConsensusTask], max_per_group: usize) -> Vec<TaskGroup> {
    if max_per_group == 0 {
        group_tasks_unlimited(tasks)
    } else {
        group_tasks_limited(tasks, max_per_group)
    }
}

fn group_tasks_unlimited(tasks: &[ConsensusTask]) -> Vec<TaskGroup> {
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

fn group_tasks_limited(tasks: &[ConsensusTask], max_per_group: usize) -> Vec<TaskGroup> {
    debug_assert!(max_per_group > 0);
    let mut group_index: HashMap<TaskGroupKey, usize> = HashMap::new();
    let mut builders: Vec<TaskGroupBuilder> = Vec::new();
    let chunk_capacity = max_per_group.min(128);

    for (idx, task) in tasks.iter().enumerate() {
        let key = task_group_key(task);
        let bidx = if let Some(&bidx) = group_index.get(&key) {
            bidx
        } else {
            let bidx = builders.len();
            group_index.insert(key.clone(), bidx);
            builders.push(TaskGroupBuilder {
                key,
                chunks: vec![Vec::with_capacity(chunk_capacity)],
            });
            bidx
        };

        let builder = &mut builders[bidx];
        if builder
            .chunks
            .last()
            .is_some_and(|chunk| chunk.len() == max_per_group)
        {
            builder.chunks.push(Vec::with_capacity(chunk_capacity));
        }
        builder
            .chunks
            .last_mut()
            .expect("limited group builder always has a chunk")
            .push(idx);
    }

    let n_groups = builders.iter().map(|builder| builder.chunks.len()).sum();
    let mut groups = Vec::with_capacity(n_groups);
    for mut builder in builders {
        let last = builder.chunks.pop();
        for indices in builder.chunks {
            groups.push(TaskGroup {
                key: builder.key.clone(),
                indices,
            });
        }
        if let Some(indices) = last {
            groups.push(TaskGroup {
                key: builder.key,
                indices,
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

fn active_samples_need_gt_fallback(
    gt_bits: &crate::vcf_store::BiallelicPhasedGtBits,
    active_words_by_hap: &[Vec<u64>; 2],
    active_word_indices_by_hap: &[Vec<usize>; 2],
) -> bool {
    let fallback_words = gt_bits.fallback_words();
    for hap_idx in 0..2 {
        let active_words = &active_words_by_hap[hap_idx];
        for &word_idx in &active_word_indices_by_hap[hap_idx] {
            if fallback_words[word_idx] & active_words[word_idx] != 0 {
                return true;
            }
        }
    }
    false
}

#[inline]
fn batch_buffer_count(output_map: &[(usize, usize)]) -> Option<usize> {
    output_map
        .iter()
        .map(|&(_, batch_idx)| batch_idx)
        .max()
        .map(|idx| idx + 1)
}

fn finish_batch_outputs(
    tasks: &[ConsensusTask],
    output_map: &[(usize, usize)],
    mut buffers: Vec<Vec<u8>>,
) -> Vec<(usize, ConsensusResult)> {
    let mut remaining = vec![0usize; buffers.len()];
    for &(_, batch_idx) in output_map {
        debug_assert!(batch_idx < buffers.len());
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
    out
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
fn byte_case_flags(byte: u8) -> u8 {
    let mut flags = 0u8;
    if byte.is_ascii_lowercase() {
        flags |= ALLELE_HAS_ASCII_LOWER;
    }
    if byte.is_ascii_uppercase() {
        flags |= ALLELE_HAS_ASCII_UPPER;
    }
    flags
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

    fn setup_phased_batch_with_unphased_inactive(
        name: &str,
    ) -> (std::path::PathBuf, std::path::PathBuf) {
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
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
             chr1\t2\t.\tC\tG\t.\t.\t.\tGT\t0|1\t1|0\t0/1\n\
             chr1\t3\t.\tG\tT\t.\t.\t.\tGT\t1|0\t0|1\t0/0\n",
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

    fn setup_phased_batch_partial_missing(name: &str) -> (std::path::PathBuf, std::path::PathBuf) {
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
             chr1\t2\t.\tC\tG\t.\t.\t.\tGT\t.|1\t0|.\n",
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
    fn biallelic_batch_flavor_specializes_option_combinations() {
        let (ref_fa, _vcf) = setup("batch_flavor");
        let new_engine = |opts: EngineOptions| {
            ConsensusEngine::new(RefIndex::load(&ref_fa).unwrap(), HashMap::new(), opts)
        };

        let engine = new_engine(EngineOptions::default());
        assert_eq!(engine.biallelic_batch_flavor(), BatchExecutionFlavor::Plain);

        let engine = new_engine(EngineOptions {
            missing: Some(b'N'),
            ..Default::default()
        });
        assert_eq!(
            engine.biallelic_batch_flavor(),
            BatchExecutionFlavor::Missing { missing: b'N' }
        );

        let engine = new_engine(EngineOptions {
            absent: Some(b'-'),
            ..Default::default()
        });
        assert_eq!(
            engine.biallelic_batch_flavor(),
            BatchExecutionFlavor::Absent { absent: b'-' }
        );

        let engine = new_engine(EngineOptions {
            absent: Some(b'-'),
            missing: Some(b'N'),
            ..Default::default()
        });
        assert_eq!(
            engine.biallelic_batch_flavor(),
            BatchExecutionFlavor::AbsentMissing {
                absent: b'-',
                missing: b'N'
            }
        );
    }

    #[test]
    fn env_disable_biallelic_batch_fastpath_returns_none_when_enabled() {
        if !biallelic_batch_fastpath_disabled() {
            return;
        }
        let (ref_fa, vcf) = setup_phased_batch("env_disable_biallelic_batch");
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

        let groups = group_tasks(&tasks, 0);
        let vcf = engine.vcfs.get("chr1").unwrap();
        let ref_seq = engine.ref_index.fetch_1based("chr1", 1, 8).unwrap();
        let (records, plan) = vcf.plan_query_set(
            "chr1",
            0,
            7,
            engine.opts.regions_overlap,
            engine.plan_options(),
        );

        assert!(engine
            .try_run_biallelic_phased_batch(&tasks, &groups[0], vcf, &ref_seq, 0, &records, &plan)
            .is_none());
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
        let groups = group_tasks(&tasks, 0);
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
        let duplicate_groups = group_tasks(&duplicate, 0);
        assert!(has_duplicate_exec_keys(
            &duplicate,
            &duplicate_groups[0].indices
        ));
    }

    #[test]
    fn group_tasks_splits_large_groups_when_limited() {
        let make_task = |idx: usize, start: i64, vcf_key: &str| ConsensusTask {
            chr: "chr1".into(),
            start,
            end: start + 7,
            vcf_key: vcf_key.into(),
            gene_id: format!("G{}", idx),
            sample: Some(format!("S{}", idx)),
            haplotype: Some("1".into()),
        };

        let tasks: Vec<_> = (0..10).map(|idx| make_task(idx, 1, "a")).collect();
        let unlimited = group_tasks(&tasks, 0);
        assert_eq!(unlimited.len(), 1);
        assert_eq!(unlimited[0].indices, (0..10).collect::<Vec<_>>());

        let groups = group_tasks(&tasks, 3);
        assert_eq!(groups.len(), 4);
        assert_eq!(groups[0].indices, vec![0, 1, 2]);
        assert_eq!(groups[1].indices, vec![3, 4, 5]);
        assert_eq!(groups[2].indices, vec![6, 7, 8]);
        assert_eq!(groups[3].indices, vec![9]);
        assert!(groups.iter().all(|group| group.indices.len() <= 3));

        let interleaved = vec![
            make_task(0, 1, "a"),
            make_task(1, 9, "a"),
            make_task(2, 1, "a"),
            make_task(3, 1, "a"),
            make_task(4, 1, "a"),
        ];
        let groups = group_tasks(&interleaved, 2);
        assert_eq!(groups.len(), 3);
        assert_eq!(groups[0].indices, vec![0, 2]);
        assert_eq!(groups[1].indices, vec![3, 4]);
        assert_eq!(groups[2].indices, vec![1]);
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
    fn consensus_stats_consume_results_without_returning_sequences() {
        let (ref_fa, vcf) = setup("stats");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(ref_fa, vcf_map, EngineOptions::default()).unwrap();
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

        let many = engine.consensus_many_stats(tasks.clone(), 4).unwrap();
        assert_eq!(many.as_tuple(), (20, 160, 8, 8));

        let iter = engine
            .consensus_iter_stats(tasks, 1, true, false, 4)
            .unwrap();
        assert_eq!(iter, many);
    }

    #[test]
    fn consensus_profile_reports_dispatch_and_compile_counters() {
        let (ref_fa, vcf) = setup_phased_batch("profile_counters");
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
                haplotype: Some("1pIu".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H2".into(),
                sample: Some("S1".into()),
                haplotype: Some("2pIu".into()),
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

        let profile = engine.consensus_many_profile(tasks, 2).unwrap();
        assert_eq!(profile.run.as_tuple(), (4, 32, 8, 8));
        assert_eq!(profile.runtime.regions_total, 1);
        assert_eq!(profile.runtime.tasks_total, 4);
        assert_eq!(profile.runtime.records_seen, 2);
        assert_eq!(profile.runtime.alloc_bytes, 32);
        assert_eq!(
            profile
                .runtime
                .lane_count(FastPathLane::BiallelicPhasedBatch),
            1
        );
        assert_eq!(profile.runtime.fallback_records, 0);
        assert!(profile
            .summary_lines()
            .contains(&"runtime.lane.BiallelicPhasedBatch=1".to_string()));
        assert!(profile
            .summary_lines()
            .iter()
            .any(|line| line.starts_with("run.seq_per_sec=")));
        assert!(profile
            .summary_lines()
            .iter()
            .any(|line| line.starts_with("runtime.records_per_sec=")));
        assert!(profile
            .summary_lines()
            .contains(&"runtime.alloc_bytes=32".to_string()));

        let compile_lines = engine.compile_stats_lines();
        assert!(compile_lines.contains(&"vcf.chr1.records_total=2".to_string()));
        assert!(compile_lines.contains(&"vcf.chr1.biallelic_gt_bitset_records=2".to_string()));
        assert!(compile_lines
            .iter()
            .any(|line| line.starts_with("vcf.chr1.allele_op.SameLen=")));
    }

    #[test]
    fn consensus_split_preserves_results_and_stats() {
        let (ref_fa, vcf) = setup_phased_batch("split_preserves_results");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let unlimited =
            ConsensusEngine::load(ref_fa.clone(), vcf_map.clone(), EngineOptions::default())
                .unwrap();
        let limited = ConsensusEngine::load(
            ref_fa,
            vcf_map,
            EngineOptions {
                max_tasks_per_group: 1,
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

        let expected = unlimited.consensus_many(tasks.clone(), 2);
        let got = limited.consensus_many(tasks.clone(), 2);
        assert_eq!(got.len(), expected.len());
        for (got, expected) in got.iter().zip(&expected) {
            assert_eq!(got.gene_id, expected.gene_id);
            assert_eq!(got.sample, expected.sample);
            assert_eq!(got.haplotype, expected.haplotype);
            assert_eq!(got.seq, expected.seq);
            assert_eq!(got.error, expected.error);
        }

        let expected_stats = unlimited.consensus_many_stats(tasks.clone(), 2).unwrap();
        assert_eq!(
            limited.consensus_many_stats(tasks.clone(), 2).unwrap(),
            expected_stats
        );
        assert_eq!(
            limited
                .consensus_iter_stats(tasks, 2, true, false, 2)
                .unwrap(),
            expected_stats
        );
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
    fn consensus_iter_handles_split_groups_ordered_and_prefetch_zero() {
        let (ref_fa, vcf) = setup_two_regions("iter_split_groups");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(
            ref_fa,
            vcf_map,
            EngineOptions {
                max_tasks_per_group: 1,
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

        let ordered: Vec<_> = engine
            .consensus_iter(tasks.clone(), 1, true, true, 2)
            .collect();
        assert_eq!(ordered.len(), 3);
        assert_eq!(ordered[0].0, 0);
        assert_eq!(ordered[1].0, 1);
        assert_eq!(ordered[2].0, 2);
        assert_eq!(ordered[0].1.seq, b"AGGTACGT");
        assert_eq!(ordered[1].1.seq, b"AAGT");
        assert_eq!(ordered[2].1.seq, b"AGGTACGT");

        let mut unordered: Vec<_> = engine.consensus_iter(tasks, 0, false, false, 2).collect();
        unordered.sort_by_key(|(idx, _)| *idx);
        assert_eq!(unordered.len(), 3);
        assert_eq!(unordered[0].0, 0);
        assert_eq!(unordered[1].0, 1);
        assert_eq!(unordered[2].0, 2);
        assert_eq!(unordered[0].1.seq, b"AGGTACGT");
        assert_eq!(unordered[1].1.seq, b"AAGT");
        assert_eq!(unordered[2].1.seq, b"AGGTACGT");
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

        let groups = group_tasks(&tasks, 0);
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
    fn biallelic_phased_batch_ignores_unphased_inactive_samples() {
        let (ref_fa, vcf) = setup_phased_batch_with_unphased_inactive("phased_batch_unphased");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(ref_fa, vcf_map, EngineOptions::default()).unwrap();
        assert_eq!(
            engine
                .vcfs
                .get("chr1")
                .unwrap()
                .compile_stats()
                .biallelic_gt_bitset_records,
            2
        );

        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H1".into(),
                sample: Some("S1".into()),
                haplotype: Some("1pIu".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S1H2".into(),
                sample: Some("S1".into()),
                haplotype: Some("2pIu".into()),
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

        let profile = engine.consensus_many_profile(tasks, 2).unwrap();
        assert_eq!(
            profile
                .runtime
                .lane_count(FastPathLane::BiallelicPhasedBatch),
            1
        );
        assert_eq!(profile.runtime.fallback_records, 0);
    }

    #[test]
    fn biallelic_phased_batch_falls_back_for_active_unphased_samples() {
        let (ref_fa, vcf) =
            setup_phased_batch_with_unphased_inactive("phased_batch_active_unphased");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(ref_fa, vcf_map, EngineOptions::default()).unwrap();
        let tasks = vec![
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S3H1".into(),
                sample: Some("S3".into()),
                haplotype: Some("1pIu".into()),
            },
            ConsensusTask {
                chr: "chr1".into(),
                start: 1,
                end: 8,
                vcf_key: "chr1".into(),
                gene_id: "S3H2".into(),
                sample: Some("S3".into()),
                haplotype: Some("2pIu".into()),
            },
        ];

        let profile = engine.consensus_many_profile(tasks.clone(), 2).unwrap();
        assert_eq!(
            profile
                .runtime
                .lane_count(FastPathLane::BiallelicPhasedBatch),
            0
        );
        assert_eq!(profile.runtime.lane_count(FastPathLane::SameLenIupac), 2);
        assert_eq!(profile.runtime.fallback_records, 0);

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"ASGTACGT");
        assert_eq!(results[1].seq, b"ASGTACGT");
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

        let groups = group_tasks(&tasks, 0);
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

        let groups = group_tasks(&tasks, 0);
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

        let groups = group_tasks(&tasks, 0);
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

        let groups = group_tasks(&tasks, 0);
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

        let groups = group_tasks(&tasks, 0);
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

        let groups = group_tasks(&tasks, 0);
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

        let groups = group_tasks(&tasks, 0);
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

        let groups = group_tasks(&tasks, 0);
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

        let groups = group_tasks(&tasks, 0);
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

        let groups = group_tasks(&tasks, 0);
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

        let groups = group_tasks(&tasks, 0);
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
    fn biallelic_phased_batch_lane_marks_missing_snv() {
        let (ref_fa, vcf) = setup_phased_batch_missing("phased_batch_missing_mark_snv");
        let mut vcf_map = HashMap::new();
        vcf_map.insert("chr1".to_string(), vcf);
        let engine = ConsensusEngine::load(
            ref_fa,
            vcf_map,
            EngineOptions {
                missing: Some(b'?'),
                mark_snv: Some(b'#'),
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
                gene_id: "S2H2".into(),
                sample: Some("S2".into()),
                haplotype: Some("2".into()),
            },
        ];

        let groups = group_tasks(&tasks, 0);
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
            .expect("batch lane should accept missing GT with -M and mark_snv");
        let mut by_idx = vec![Vec::new(); batch.len()];
        for (idx, result) in batch {
            by_idx[idx] = result.seq;
        }
        assert_eq!(by_idx[0], b"A#GTACGT");
        assert_eq!(by_idx[1], b"A#GTACGT");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"A#GTACGT");
        assert_eq!(results[1].seq, b"A#GTACGT");
    }

    #[test]
    fn biallelic_phased_batch_lane_handles_partial_missing_per_haplotype() {
        let (ref_fa, vcf) = setup_phased_batch_partial_missing("phased_batch_partial_missing");
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

        let groups = group_tasks(&tasks, 0);
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
            .expect("batch lane should accept phased partial missing GT with -M");
        let mut by_idx = vec![Vec::new(); batch.len()];
        for (idx, result) in batch {
            by_idx[idx] = result.seq;
        }
        assert_eq!(by_idx[0], b"A?GTACGT");
        assert_eq!(by_idx[1], b"AGGTACGT");
        assert_eq!(by_idx[2], b"ACGTACGT");
        assert_eq!(by_idx[3], b"A?GTACGT");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"A?GTACGT");
        assert_eq!(results[1].seq, b"AGGTACGT");
        assert_eq!(results[2].seq, b"ACGTACGT");
        assert_eq!(results[3].seq, b"A?GTACGT");
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

        let groups = group_tasks(&tasks, 0);
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
        assert_eq!(by_idx[0], b"A?GTACGT");
        assert_eq!(by_idx[1], b"A?GTACGT");
        assert_eq!(by_idx[2], b"ACGTACGT");
        assert_eq!(by_idx[3], b"ATTTACGT");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"A?GTACGT");
        assert_eq!(results[1].seq, b"A?GTACGT");
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

        let groups = group_tasks(&tasks, 0);
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
        assert_eq!(by_idx[0], b"N?NNNNNN");
        assert_eq!(by_idx[1], b"NCGNNNNN");
        assert_eq!(by_idx[2], b"NTTNNNNN");

        let results = engine.consensus_many(tasks, 2);
        assert_eq!(results[0].seq, b"N?NNNNNN");
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

        let groups = group_tasks(&tasks, 0);
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
