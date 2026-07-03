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

use crate::engine::{ConsensusEngine, ConsensusTask, EngineOptions};
use crate::mask::MaskWith;
use pyo3::exceptions::{PyIOError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

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
    seq: Vec<u8>,
    #[pyo3(get)]
    chain: Option<String>,
}

#[pyclass(module = "pyconsensus._engine", name = "_ConsensusEngine", subclassable)]
struct PyConsensusEngine {
    inner: Arc<ConsensusEngine>,
}

#[pymethods]
impl PyConsensusEngine {
    /// Construct + eagerly load ref and all VCFs.
    #[new]
    #[pyo3(signature = (ref_path, vcfs, iupac_codes=false, missing=None, absent=None, mark_del=None, mark_ins=None, mark_snv=None, mask=None, mask_with="N".to_string(), chain=false, regions_overlap=1u8))]
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
        };
        let vcf_paths: HashMap<String, PathBuf> = vcfs
            .into_iter()
            .map(|(k, v)| (k, PathBuf::from(v)))
            .collect();
        let engine =
            ConsensusEngine::load(ref_path, vcf_paths, opts).map_err(|e| PyIOError::new_err(e))?;
        Ok(PyConsensusEngine {
            inner: Arc::new(engine),
        })
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
        let mut items: Vec<PyConsensusResult> = Vec::with_capacity(results.len());
        for r in results {
            if let Some(e) = r.error {
                return Err(PyRuntimeError::new_err(format!("{}: {}", r.gene_id, e)));
            }
            items.push(PyConsensusResult {
                gene_id: r.gene_id,
                sample: r.sample,
                haplotype: r.haplotype,
                seq: r.seq,
                chain: r.chain,
            });
        }
        Ok(items.into_pyobject(py)?.into_any())
    }

    /// Launch a lazy iterator over `Task` objects using a producer-consumer
    /// model.
    ///
    /// * `prefetch_steps` — 0 = no prefetch (next submits 1, blocks on it);
    ///   N>0 = keep N tasks in flight.
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
    #[pyo3(signature = (tasks, prefetch_steps=0, warmup=false, ordered=false, threads=1))]
    fn consensus_iter(
        &self,
        tasks: Vec<PyRef<'_, PyTask>>,
        prefetch_steps: usize,
        warmup: bool,
        ordered: bool,
        threads: usize,
    ) -> PyConsensusIter {
        let engine = self.inner.clone();
        let tasks: Vec<ConsensusTask> = tasks.iter().map(|t| t.to_inner()).collect();
        let iter = engine.consensus_iter(tasks, prefetch_steps, warmup, ordered, threads);
        PyConsensusIter { inner: Some(iter) }
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
                if let Some(e) = r.error {
                    return Err(PyRuntimeError::new_err(format!("{}: {}", r.gene_id, e)));
                }
                let obj = PyConsensusResult {
                    gene_id: r.gene_id,
                    sample: r.sample,
                    haplotype: r.haplotype,
                    seq: r.seq,
                    chain: r.chain,
                };
                let result_obj = obj.into_pyobject(py)?.into_any().unbind();
                let idx_obj = idx.into_pyobject(py)?.into_any().unbind();
                let tup = PyTuple::new(py, &[idx_obj, result_obj])?
                    .into_any()
                    .unbind();
                Ok(Some(tup))
            }
        }
    }
}

#[pymodule]
fn _engine(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyTask>()?;
    m.add_class::<PyConsensusResult>()?;
    m.add_class::<PyConsensusEngine>()?;
    m.add_class::<PyConsensusIter>()?;
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
