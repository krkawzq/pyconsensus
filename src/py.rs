//! py — PyO3 bindings (feature `python`).
//!
//! Exposes the *private* `_ConsensusEngine` / `_ConsensusIter` Rust classes plus
//! public `Task` and `ConsensusResult` dataclass-style objects. The public
//! `ConsensusEngine` wrapper and the `build_tasks` helper live in the
//! pure-Python `pyconsensus.engine` module, which layers Python-side logic on
//! top of this private extension module.
//!
//! The iterator's `__next__` releases the GIL while blocking on the result
//! channel, so Rust worker threads can produce in parallel with Python
//! consumption.
//!
//! Each produced item is a `ConsensusResult` with named fields (`gene_id`,
//! `sample`, `haplotype`, `seq`, `chain`) — `chain` is always present (None
//! when chain output is disabled). `consensus_many` returns them in input
//! order as a flat list (list index == input position). The lazy iterator
//! yields `(idx, ConsensusResult)` tuples — `idx` is the input position, needed
//! to re-pair results in unordered completion mode.

#![cfg(feature = "python")]

use crate::apply::ApplyErrorMode;
use crate::engine::{
    cache_dedupe_key, compile_thread_count, thread_pool, ConsensusEngine, ConsensusResult,
    ConsensusTask, EngineOptions,
};
use crate::logging::{
    ensure_default_htslib_log_level, htslib_log_level, set_htslib_log_level, HtsLogLevel,
    LogControl, LogLevel,
};
use crate::mask::MaskWith;
use crate::vcf_store::VcfStore;
use pyo3::exceptions::{PyIOError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyList, PyTuple};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

const APPLY_WARNING_LIMIT: usize = 20;
static APPLY_WARNING_COUNT: AtomicUsize = AtomicUsize::new(0);

fn warn_apply_failure(log: &LogControl, context: &str, gene_id: &str, error: &str) {
    if !log.enabled(LogLevel::Warn) {
        return;
    }
    let n = APPLY_WARNING_COUNT.fetch_add(1, Ordering::Relaxed);
    if n < APPLY_WARNING_LIMIT {
        eprintln!("pyconsensus warning: {context}: {gene_id}: {error}");
    } else if n == APPLY_WARNING_LIMIT {
        eprintln!(
            "pyconsensus warning: further apply failures suppressed after {APPLY_WARNING_LIMIT} messages"
        );
    }
}

fn warn_apply_failure_summary(
    log: &LogControl,
    context: &str,
    failed: usize,
    first_error: Option<&str>,
) {
    if failed == 0 {
        return;
    }
    let detail = first_error.unwrap_or("no detail");
    warn_apply_failure(
        log,
        context,
        "summary",
        &format!("{failed} failed task(s); first: {detail}"),
    );
}

fn py_consensus_result(r: ConsensusResult) -> PyConsensusResult {
    let failed = r.error.is_some();
    PyConsensusResult {
        gene_id: r.gene_id,
        sample: r.sample,
        haplotype: r.haplotype,
        seq: if failed { None } else { Some(r.seq) },
        chain: if failed { None } else { r.chain },
        error: r.error,
    }
}

fn convert_result(
    context: &str,
    mode: ApplyErrorMode,
    log: &LogControl,
    r: ConsensusResult,
) -> PyResult<PyConsensusResult> {
    if let Some(error) = r.error.as_deref() {
        if mode.is_error() {
            return Err(PyRuntimeError::new_err(format!("{}: {}", r.gene_id, error)));
        }
        warn_apply_failure(log, context, &r.gene_id, error);
    }
    Ok(py_consensus_result(r))
}

/// Parse a cli-style mask-with string ("uc"/"lc"/single char).
fn parse_mask_with(s: &str) -> PyResult<MaskWith> {
    let lower = s.to_ascii_lowercase();
    match lower.as_str() {
        "uc" => Ok(MaskWith::Uc),
        "lc" => Ok(MaskWith::Lc),
        _ if s.chars().count() == 1 => Ok(MaskWith::Char(s.as_bytes()[0])),
        _ => Err(PyValueError::new_err(format!(
            "mask_with must be 'uc', 'lc', or a single char, got {:?}",
            s
        ))),
    }
}

fn parse_log_level(s: &str) -> PyResult<LogLevel> {
    LogLevel::parse(s).ok_or_else(|| {
        PyValueError::new_err(format!(
            "log_level must be one of 'off', 'error', 'warn', 'info', or 'debug', got {:?}",
            s
        ))
    })
}

fn parse_htslib_log_level(s: &str) -> PyResult<HtsLogLevel> {
    HtsLogLevel::parse(s).ok_or_else(|| {
        PyValueError::new_err(format!(
            "htslib log_level must be one of 'off', 'error', 'warn', 'info', 'debug', or 'trace', got {:?}",
            s
        ))
    })
}

#[pyfunction]
fn get_htslib_log_level() -> String {
    htslib_log_level().as_str().to_string()
}

#[pyfunction(name = "set_htslib_log_level")]
fn set_htslib_log_level_py(level: &str) -> PyResult<()> {
    set_htslib_log_level(parse_htslib_log_level(level)?);
    Ok(())
}

/// One consensus production request, exposed to Python as a dataclass-style
/// object with get/set fields. Constructed directly from Python
/// (`Task("chr1", 1, 8, "chr1", "G1", sample=..., haplotype=...)`); passed
/// verbatim to `consensus_many` / `consensus_iter`.
#[pyclass(module = "pyconsensus._engine", name = "Task")]
#[derive(Clone)]
struct PyTask {
    #[pyo3(get, set)]
    chr: String,
    #[pyo3(get, set)]
    start: i64,
    #[pyo3(get, set)]
    end: i64,
    #[pyo3(get, set)]
    vcf_key: String,
    #[pyo3(get, set)]
    gene_id: String,
    #[pyo3(get, set)]
    sample: Option<String>,
    #[pyo3(get, set)]
    haplotype: Option<String>,
}

impl PyTask {
    /// Convert to the internal Rust task type used by the engine core.
    fn to_inner(&self) -> ConsensusTask {
        ConsensusTask {
            chr: self.chr.clone(),
            start: self.start,
            end: self.end,
            vcf_key: self.vcf_key.clone(),
            gene_id: self.gene_id.clone(),
            sample: self.sample.clone(),
            haplotype: self.haplotype.clone(),
        }
    }
}

#[pymethods]
impl PyTask {
    /// `Task(chr, start, end, vcf_key, gene_id, sample=None, haplotype=None)`.
    ///
    /// `sample` / `haplotype` optional:
    /// * both None — apply all ALT (or IUPAC from REF/ALT when iupac_codes=True);
    /// * sample set, haplotype None — IUPAC across that sample's GTs;
    /// * both set — single-sample haplotype selection (e.g. "1pIu" / "2pIu").
    #[new]
    #[pyo3(signature = (chr, start, end, vcf_key, gene_id, sample=None, haplotype=None))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        chr: String,
        start: i64,
        end: i64,
        vcf_key: String,
        gene_id: String,
        sample: Option<String>,
        haplotype: Option<String>,
    ) -> Self {
        PyTask {
            chr,
            start,
            end,
            vcf_key,
            gene_id,
            sample,
            haplotype,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Task(chr={:?}, start={}, end={}, vcf_key={:?}, gene_id={:?}, sample={:?}, haplotype={:?})",
            self.chr, self.start, self.end, self.vcf_key, self.gene_id, self.sample, self.haplotype
        )
    }
}

/// One produced consensus result, exposed to Python with named fields. `idx`
/// is NOT a field here — it is the input position, returned alongside the
/// result by the lazy iterator (`(idx, ConsensusResult)`). `consensus_many`
/// returns results in input order, so the list index is the position.
#[pyclass(module = "pyconsensus._engine", name = "ConsensusResult")]
struct PyConsensusResult {
    #[pyo3(get)]
    gene_id: String,
    #[pyo3(get)]
    sample: Option<String>,
    #[pyo3(get)]
    haplotype: Option<String>,
    #[pyo3(get)]
    seq: Option<Vec<u8>>,
    #[pyo3(get)]
    chain: Option<String>,
    #[pyo3(get)]
    error: Option<String>,
}

/// One `.cvcf` cache build result (dataclass-style, Rust-backed).
///
/// Returned by `ConsensusEngine.build_cache`. `status` is one of:
/// `"hit"` (valid cache read as-is), `"built"` (no prior cache, parsed +
/// written), `"rebuilt"` (prior cache failed validation, reparsed +
/// rewritten), `"forced"` (`force=True`, cache ignored, reparsed + rewritten).
#[pyclass(module = "pyconsensus._engine", name = "CacheResult")]
struct PyCacheResult {
    #[pyo3(get)]
    path: String,
    #[pyo3(get)]
    cache_path: String,
    #[pyo3(get)]
    status: String,
    #[pyo3(get)]
    records: usize,
    #[pyo3(get)]
    samples: usize,
    #[pyo3(get)]
    cache_mb: f64,
    #[pyo3(get)]
    elapsed_sec: f64,
}

#[pyclass(module = "pyconsensus._engine", name = "_ConsensusEngine", subclass)]
struct PyConsensusEngine {
    inner: Arc<ConsensusEngine>,
}

#[pymethods]
impl PyConsensusEngine {
    /// Construct + eagerly load ref and all VCFs.
    #[new]
    #[pyo3(signature = (ref_path, vcfs, iupac_codes=false, missing=None, absent=None, mark_del=None, mark_ins=None, mark_snv=None, mask=None, mask_with="N".to_string(), chain=false, regions_overlap=0u8, max_tasks_per_group=0usize, compile_threads=None, log_level="info".to_string()))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        ref_path: String,
        vcfs: HashMap<String, String>,
        iupac_codes: bool,
        missing: Option<String>,
        absent: Option<String>,
        mark_del: Option<String>,
        mark_ins: Option<String>,
        mark_snv: Option<String>,
        mask: Option<String>,
        mask_with: String,
        chain: bool,
        regions_overlap: u8,
        max_tasks_per_group: usize,
        compile_threads: Option<usize>,
        log_level: String,
    ) -> PyResult<Self> {
        let one_char = |s: &Option<String>, name: &str| -> PyResult<Option<u8>> {
            match s {
                None => Ok(None),
                Some(s) if s.chars().count() == 1 => Ok(Some(s.as_bytes()[0])),
                _ => Err(PyValueError::new_err(format!(
                    "{} must be a single char",
                    name
                ))),
            }
        };
        // mark_ins / mark_snv accept "uc"/"lc"/char
        let mark_val = |s: &Option<String>, name: &str| -> PyResult<Option<u8>> {
            match s {
                None => Ok(None),
                Some(v) => match v.as_str() {
                    "uc" => Ok(Some(1)),
                    "lc" => Ok(Some(2)),
                    other if other.chars().count() == 1 => Ok(Some(other.as_bytes()[0])),
                    _ => Err(PyValueError::new_err(format!(
                        "{} must be 'uc','lc',or a char",
                        name
                    ))),
                },
            }
        };

        let opts = EngineOptions {
            iupac_codes,
            missing: one_char(&missing, "missing")?,
            absent: one_char(&absent, "absent")?,
            mark_del: one_char(&mark_del, "mark_del")?,
            mark_ins: mark_val(&mark_ins, "mark_ins")?,
            mark_snv: mark_val(&mark_snv, "mark_snv")?,
            mask: mask.as_ref().map(PathBuf::from),
            mask_with: parse_mask_with(&mask_with)?,
            chain,
            regions_overlap,
            apply_error_mode: ApplyErrorMode::from_env(),
            compile_threads,
            log_level: parse_log_level(&log_level)?,
            max_tasks_per_group,
        };
        let vcf_paths: HashMap<String, PathBuf> = vcfs
            .into_iter()
            .map(|(k, v)| (k, PathBuf::from(v)))
            .collect();
        let engine =
            ConsensusEngine::load(ref_path, vcf_paths, opts).map_err(PyIOError::new_err)?;
        Ok(PyConsensusEngine {
            inner: Arc::new(engine),
        })
    }

    /// Pre-build `.cvcf` caches for a list of VCF files without loading the
    /// reference or constructing the engine.
    ///
    /// `paths` — list of VCF/BCF file paths (plain or bgzipped). Each is
    /// loaded with the same cache logic as the constructor: an existing valid
    /// cache is read as-is (`status="hit"`); a missing or invalid cache is
    /// reparsed and rewritten (`status="built"` / `"rebuilt"`). Paths that
    /// resolve to the same cache file (after canonicalization) are loaded
    /// only once; duplicates are skipped.
    ///
    /// `compile_threads` — size of the rayon pool that loads VCFs in parallel
    /// (one thread per unique VCF; a single VCF is parsed single-threaded).
    /// `None` uses available parallelism capped at the unique VCF count,
    /// matching the constructor's `compile_threads`.
    ///
    /// `force` — if true, ignore any existing cache and rebuild unconditionally
    /// (`status="forced"`). When false, an invalid cache is still rebuilt.
    ///
    /// Returns one `CacheResult` per unique input VCF, in first-seen order.
    /// On any VCF failure, raises `IOError` with the offending path.
    #[staticmethod]
    #[pyo3(signature = (paths, compile_threads=None, force=false, log_level="info".to_string()))]
    fn build_cache(
        paths: Vec<String>,
        compile_threads: Option<usize>,
        force: bool,
        log_level: String,
    ) -> PyResult<Vec<PyCacheResult>> {
        let log = LogControl::new(parse_log_level(&log_level)?);

        // Deduplicate by canonical cache path (same logic as the constructor's
        // `group_vcf_loads`): paths resolving to the same `.cvcf` are loaded
        // only once. Keep first-seen order for a stable result list.
        let mut seen: HashSet<PathBuf> = HashSet::new();
        let mut unique: Vec<PathBuf> = Vec::new();
        for p in paths {
            let path = PathBuf::from(p);
            let dedup_key = cache_dedupe_key(&path);
            if seen.insert(dedup_key) {
                unique.push(path);
            } else if log.enabled(LogLevel::Warn) {
                eprintln!(
                    "pyconsensus warning: build_cache_skip_duplicate path={}",
                    path.display()
                );
            }
        }

        let load_threads = compile_thread_count(compile_threads, unique.len());
        if log.enabled(LogLevel::Info) {
            eprintln!(
                "pyconsensus info: build_cache_start vcf_count={} compile_threads={} force={}",
                unique.len(),
                load_threads,
                force
            );
        }

        if unique.is_empty() {
            return Ok(Vec::new());
        }

        // Load in parallel; carry the first-seen index so we can restore order
        // after par_iter (which completes in arbitrary order).
        let pool = thread_pool(load_threads);
        let loaded: Result<Vec<(usize, PyCacheResult)>, String> = pool.install(|| {
            unique
                .into_par_iter()
                .enumerate()
                .map(
                    |(idx, path)| match VcfStore::load_with_outcome(&path, None, &log, force) {
                        Ok((store, outcome)) => Ok((
                            idx,
                            PyCacheResult {
                                path: store.path().to_string_lossy().into_owned(),
                                cache_path: store.cache_path().to_string_lossy().into_owned(),
                                status: outcome.status.to_string(),
                                records: outcome.records,
                                samples: outcome.samples,
                                cache_mb: outcome.cache_mb,
                                elapsed_sec: outcome.elapsed_sec,
                            },
                        )),
                        Err(err) => Err(format!("VCF '{}': {}", path.display(), err)),
                    },
                )
                .collect()
        });
        let mut loaded = loaded.map_err(PyIOError::new_err)?;
        loaded.sort_by_key(|(idx, _)| *idx);
        Ok(loaded.into_iter().map(|(_, r)| r).collect())
    }

    #[getter]
    fn log_level(&self) -> String {
        self.inner.log_level().as_str().to_string()
    }

    #[setter]
    fn set_log_level(&self, level: &str) -> PyResult<()> {
        self.inner.set_log_level(parse_log_level(level)?);
        Ok(())
    }

    /// Run a flat list of `Task` objects in parallel, returning
    /// `ConsensusResult` objects in input order (list index == input position).
    ///
    /// `threads` — size of the rayon pool created for this call. A fresh pool
    /// is built per call; the engine itself holds no compute resources, only
    /// the preprocessed ref + VCFs.
    #[pyo3(signature = (tasks, threads=1))]
    fn consensus_many<'py>(
        &self,
        py: Python<'py>,
        tasks: Vec<PyRef<'_, PyTask>>,
        threads: usize,
    ) -> PyResult<Bound<'py, PyAny>> {
        let tasks: Vec<ConsensusTask> = tasks.iter().map(|t| t.to_inner()).collect();
        let engine = self.inner.clone();
        let results = py.allow_threads(move || engine.consensus_many(tasks, threads));
        let mode = self.inner.apply_error_mode();
        let log = self.inner.log_control();
        let mut items: Vec<PyConsensusResult> = Vec::with_capacity(results.len());
        for r in results {
            items.push(convert_result("consensus_many", mode, log.as_ref(), r)?);
        }
        Ok(items.into_pyobject(py)?.into_any())
    }

    /// Run a flat list of tasks and consume the sequences inside Rust.
    ///
    /// Returns `(n, total_len, min_len, max_len)`. This avoids building a large
    /// Python list of `ConsensusResult` objects and is intended for throughput
    /// benchmarks or sinks that discard sequence bytes.
    #[pyo3(signature = (tasks, threads=1))]
    fn consensus_many_stats(
        &self,
        py: Python<'_>,
        tasks: Vec<PyRef<'_, PyTask>>,
        threads: usize,
    ) -> PyResult<(usize, u64, usize, usize)> {
        let tasks: Vec<ConsensusTask> = tasks.iter().map(|t| t.to_inner()).collect();
        let engine = self.inner.clone();
        let stats = py
            .allow_threads(move || engine.consensus_many_stats(tasks, threads))
            .map_err(PyRuntimeError::new_err)?;
        if !self.inner.apply_error_mode().is_error() {
            let log = self.inner.log_control();
            warn_apply_failure_summary(
                log.as_ref(),
                "consensus_many_stats",
                stats.failed,
                stats.first_error.as_deref(),
            );
        }
        Ok(stats.as_tuple())
    }

    /// Run a flat list of tasks and return key-value profile lines.
    ///
    /// The profile includes throughput counters plus runtime lane/fallback
    /// counters. Sequence bytes are consumed inside Rust and are not returned.
    #[pyo3(signature = (tasks, threads=1))]
    fn consensus_many_profile(
        &self,
        py: Python<'_>,
        tasks: Vec<PyRef<'_, PyTask>>,
        threads: usize,
    ) -> PyResult<Vec<String>> {
        let tasks: Vec<ConsensusTask> = tasks.iter().map(|t| t.to_inner()).collect();
        let engine = self.inner.clone();
        let profile = py
            .allow_threads(move || engine.consensus_many_profile(tasks, threads))
            .map_err(PyRuntimeError::new_err)?;
        if !self.inner.apply_error_mode().is_error() {
            let log = self.inner.log_control();
            warn_apply_failure_summary(
                log.as_ref(),
                "consensus_many_profile",
                profile.run.failed,
                profile.run.first_error.as_deref(),
            );
        }
        Ok(profile.summary_lines())
    }

    /// Return VCF compile-time counters as key-value lines.
    fn compile_stats(&self) -> Vec<String> {
        self.inner.compile_stats_lines()
    }

    /// Launch a lazy iterator over `Task` objects using a producer-consumer
    /// model.
    ///
    /// * `prefetch_steps` — None = use `threads`; 0 = no prefetch (next
    ///   submits 1, blocks on it); N>0 = keep N region groups in flight.
    /// * `warmup` — true = start prefetching at construction; false = defer to
    ///   first `next()`.
    /// * `ordered` — true = yield in input order; false = yield in completion
    ///   order, each result carrying its input `idx`.
    /// * `threads` — rayon worker pool size for this call (fresh pool per
    ///   iterator; the engine holds no pool of its own).
    ///
    /// Each yielded item is a `(idx, ConsensusResult)` tuple — `idx` is the
    /// task's input position, needed to re-pair results when ``ordered=False``
    /// (completion order). With ``ordered=True`` idx arrives in ascending
    /// order anyway.
    #[pyo3(signature = (tasks, prefetch_steps=None, warmup=false, ordered=false, threads=1))]
    fn consensus_iter(
        &self,
        tasks: Vec<PyRef<'_, PyTask>>,
        prefetch_steps: Option<usize>,
        warmup: bool,
        ordered: bool,
        threads: usize,
    ) -> PyConsensusIter {
        let engine = self.inner.clone();
        let tasks: Vec<ConsensusTask> = tasks.iter().map(|t| t.to_inner()).collect();
        let prefetch_steps = prefetch_steps.unwrap_or_else(|| threads.max(1));
        let iter = engine.consensus_iter(tasks, prefetch_steps, warmup, ordered, threads);
        PyConsensusIter {
            inner: Some(iter),
            error_mode: self.inner.apply_error_mode(),
            log: self.inner.log_control(),
        }
    }

    /// Drive `consensus_iter` to completion and consume sequences inside Rust.
    ///
    /// Returns `(n, total_len, min_len, max_len)`.
    #[pyo3(signature = (tasks, prefetch_steps=None, warmup=false, ordered=false, threads=1))]
    fn consensus_iter_stats(
        &self,
        py: Python<'_>,
        tasks: Vec<PyRef<'_, PyTask>>,
        prefetch_steps: Option<usize>,
        warmup: bool,
        ordered: bool,
        threads: usize,
    ) -> PyResult<(usize, u64, usize, usize)> {
        let engine = self.inner.clone();
        let tasks: Vec<ConsensusTask> = tasks.iter().map(|t| t.to_inner()).collect();
        let prefetch_steps = prefetch_steps.unwrap_or_else(|| threads.max(1));
        let stats = py
            .allow_threads(move || {
                engine.consensus_iter_stats(tasks, prefetch_steps, warmup, ordered, threads)
            })
            .map_err(PyRuntimeError::new_err)?;
        if !self.inner.apply_error_mode().is_error() {
            let log = self.inner.log_control();
            warn_apply_failure_summary(
                log.as_ref(),
                "consensus_iter_stats",
                stats.failed,
                stats.first_error.as_deref(),
            );
        }
        Ok(stats.as_tuple())
    }
}

/// Lazy Python iterator over consensus results.
///
/// `__next__` releases the GIL while blocking on the completion queue, so Rust
/// worker threads keep producing in parallel with Python consumption. Each
/// call yields `(idx, ConsensusResult)`.
#[pyclass(module = "pyconsensus._engine", name = "_ConsensusIter")]
struct PyConsensusIter {
    inner: Option<crate::engine::ConsensusIter>,
    error_mode: ApplyErrorMode,
    log: Arc<LogControl>,
}

impl PyConsensusIter {
    fn inner_mode(&self) -> ApplyErrorMode {
        self.error_mode
    }
}

#[pymethods]
impl PyConsensusIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<PyObject>> {
        let inner = match self.inner.as_mut() {
            Some(i) => i,
            None => return Ok(None),
        };
        // Release the GIL while blocking on the completion queue + submitting
        // the next prefetch task (all done in Rust).
        let result = py.allow_threads(|| inner.next_blocking());
        match result {
            None => {
                self.inner = None;
                Ok(None)
            }
            Some((idx, r)) => {
                let mode = self.inner_mode();
                let obj = convert_result("consensus_iter", mode, self.log.as_ref(), r)?;
                let result_obj = obj.into_pyobject(py)?.into_any().unbind();
                let idx_obj = idx.into_pyobject(py)?.into_any().unbind();
                let tup = PyTuple::new(py, &[idx_obj, result_obj])?
                    .into_any()
                    .unbind();
                Ok(Some(tup))
            }
        }
    }

    /// Return up to `batch_size` iterator items as a Python list.
    ///
    /// This amortizes Python boundary crossings for large result streams. It
    /// returns `None` after the iterator is exhausted, matching `__next__`.
    #[pyo3(signature = (batch_size))]
    fn next_batch(&mut self, py: Python<'_>, batch_size: usize) -> PyResult<Option<PyObject>> {
        if batch_size == 0 {
            return Ok(Some(PyList::empty(py).into_any().unbind()));
        }
        let Some(inner) = self.inner.as_mut() else {
            return Ok(None);
        };
        let items = py.allow_threads(|| {
            let mut items = Vec::with_capacity(batch_size.min(1024));
            for _ in 0..batch_size {
                match inner.next_blocking() {
                    Some(item) => items.push(item),
                    None => break,
                }
            }
            items
        });
        if items.is_empty() {
            self.inner = None;
            return Ok(None);
        }
        let mut py_items = Vec::with_capacity(items.len());
        let mode = self.inner_mode();
        for (idx, r) in items {
            let obj = convert_result("consensus_iter.next_batch", mode, self.log.as_ref(), r)?;
            let result_obj = obj.into_pyobject(py)?.into_any().unbind();
            let idx_obj = idx.into_pyobject(py)?.into_any().unbind();
            py_items.push(
                PyTuple::new(py, &[idx_obj, result_obj])?
                    .into_any()
                    .unbind(),
            );
        }
        Ok(Some(PyList::new(py, py_items)?.into_any().unbind()))
    }

    /// Return up to `batch_size` `(idx, seq)` tuples as a Python list.
    ///
    /// This is the lowest-overhead Python streaming path when the caller keeps
    /// task metadata outside the result object and only needs sequence bytes.
    #[pyo3(signature = (batch_size))]
    fn next_batch_bytes(
        &mut self,
        py: Python<'_>,
        batch_size: usize,
    ) -> PyResult<Option<PyObject>> {
        if batch_size == 0 {
            return Ok(Some(PyList::empty(py).into_any().unbind()));
        }
        let Some(inner) = self.inner.as_mut() else {
            return Ok(None);
        };
        let items = py.allow_threads(|| {
            let mut items = Vec::with_capacity(batch_size.min(1024));
            for _ in 0..batch_size {
                match inner.next_blocking() {
                    Some(item) => items.push(item),
                    None => break,
                }
            }
            items
        });
        if items.is_empty() {
            self.inner = None;
            return Ok(None);
        }
        let mut py_items = Vec::with_capacity(items.len());
        let mode = self.inner_mode();
        for (idx, r) in items {
            let error = r.error.clone();
            if let Some(e) = error.as_deref() {
                if mode.is_error() {
                    return Err(PyRuntimeError::new_err(format!("{}: {}", r.gene_id, e)));
                }
                warn_apply_failure(
                    self.log.as_ref(),
                    "consensus_iter.next_batch_bytes",
                    &r.gene_id,
                    e,
                );
            }
            let idx_obj = idx.into_pyobject(py)?.into_any().unbind();
            let seq_obj = if error.is_some() {
                py.None()
            } else {
                PyBytes::new(py, &r.seq).into_any().unbind()
            };
            py_items.push(PyTuple::new(py, &[idx_obj, seq_obj])?.into_any().unbind());
        }
        Ok(Some(PyList::new(py, py_items)?.into_any().unbind()))
    }
}

#[pymodule]
fn _engine(m: &Bound<'_, PyModule>) -> PyResult<()> {
    ensure_default_htslib_log_level();
    m.add_class::<PyTask>()?;
    m.add_class::<PyConsensusResult>()?;
    m.add_class::<PyCacheResult>()?;
    m.add_class::<PyConsensusEngine>()?;
    m.add_class::<PyConsensusIter>()?;
    m.add_function(wrap_pyfunction!(get_htslib_log_level, m)?)?;
    m.add_function(wrap_pyfunction!(set_htslib_log_level_py, m)?)?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
