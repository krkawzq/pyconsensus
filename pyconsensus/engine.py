"""Public Python API for the pyconsensus engine.

A thin pure-Python layer over the private Rust extension module
`pyconsensus._engine` (which exposes `Task`, `ConsensusResult`,
`_ConsensusEngine`, `_ConsensusIter`). Python-side logic — input expansion,
validation, future pre/post-processing — lives here so it can evolve without
recompiling Rust.

Layering:
  * `Task` / `ConsensusResult` — Rust dataclass-style objects (re-exported).
  * `build_tasks(...)`          — expand regions × samples × haplotypes into a
                                  flat `list[Task]` (eager, Python-side).
  * `ConsensusEngine`           — public subclass of `_ConsensusEngine`; adds
                                  `consensus_regions` on top of the inherited
                                  `consensus_many` / `consensus_iter`.

The engine holds only preprocessed material (ref + VCFs); compute resources
(thread pool, thread count) are passed per `consensus_*` call. The Rust side
never expands a cartesian product — it receives a flat task list verbatim.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from ._engine import (
    CacheResult,
    ConsensusResult,
    Task,
    _ConsensusEngine,
    _ConsensusIter,
    __version__,
    get_htslib_log_level,
    set_htslib_log_level,
)

__all__ = [
    "CacheResult",
    "ConsensusEngine",
    "ConsensusResult",
    "Task",
    "build_cache",
    "build_tasks",
    "get_htslib_log_level",
    "set_htslib_log_level",
    "__version__",
]


def build_tasks(
    regions: Sequence[Task],
    samples: Sequence[str] | None = None,
    haplotypes: Sequence[str] | None = None,
) -> list[Task]:
    """Expand ``regions × samples × haplotypes`` into a flat task list.

    Each entry of `regions` is a `Task` whose `sample` / `haplotype` are
    ignored (treated as a region template: chr/start/end/vcf_key/gene_id). A
    `None` (or empty) dimension collapses — it yields one task with that field
    set to `None` rather than zero tasks:

    * ``build_tasks(regions)``                      -> one task/region, no sample/hap
    * ``build_tasks(regions, samples)``             -> one task/(region,sample), no hap
    * ``build_tasks(regions, samples, haplotypes)`` -> full product

    Order is region-major (region, then sample, then haplotype): all haplotypes
    of a sample for a gene are contiguous, then the next sample, then the next
    gene — matching the layout ``make_consensus_enformer_new.py`` wrote.
    """
    samp_dim: list[str | None] = list(samples) if samples else [None]
    hap_dim: list[str | None] = list(haplotypes) if haplotypes else [None]

    tasks: list[Task] = []
    for reg in regions:
        for sample in samp_dim:
            for hap in hap_dim:
                tasks.append(
                    Task(
                        reg.chr,
                        reg.start,
                        reg.end,
                        reg.vcf_key,
                        reg.gene_id,
                        sample,
                        hap,
                    )
                )
    return tasks


def build_cache(
    paths: Sequence[str],
    *,
    compile_threads: int | None = None,
    force: bool = False,
    log_level: str = "info",
) -> list[CacheResult]:
    """Pre-build ``.cvcf`` caches for a list of VCF files without loading
    the reference or constructing the engine.

    Thin wrapper over the Rust ``_ConsensusEngine.build_cache`` staticmethod.
    Each path in `paths` is a VCF/BCF file (plain or bgzipped), loaded with
    the same cache logic as the constructor: an existing valid cache is read
    as-is (``status="hit"``); a missing or invalid cache is reparsed and
    rewritten (``status="built"`` / ``"rebuilt"``). Paths resolving to the
    same ``.cvcf`` (after canonicalization) are loaded only once; duplicates
    are skipped.

    * `compile_threads` — rayon pool size for loading VCFs in parallel (one
      thread per unique VCF; a single VCF parses single-threaded). ``None``
      uses available parallelism capped at the unique VCF count, matching the
      constructor.
    * `force` — ignore any existing cache and rebuild unconditionally
      (``status="forced"``). When ``False``, an invalid cache is still rebuilt.

    Returns one :class:`CacheResult` per unique input VCF, in first-seen
    order. On any VCF failure, raises :class:`OSError` (``IOError``) naming
    the offending path.
    """
    return _ConsensusEngine.build_cache(
        list(paths),
        compile_threads=compile_threads,
        force=force,
        log_level=log_level,
    )


class ConsensusEngine(_ConsensusEngine):
    """Public engine, subclassing the private Rust `_ConsensusEngine`.

    `consensus_many` / `consensus_iter` are inherited verbatim from the Rust
    parent — no forwarding boilerplate. This subclass only adds the
    `consensus_regions` convenience method, which expands regions × samples ×
    haplotypes in Python and feeds the flat task list to the inherited
    `consensus_iter`.

    Both paths feed the Rust engine a flat task list verbatim; the engine does
    no product expansion of its own. `consensus_iter` / `consensus_regions`
    return the Rust `_ConsensusIter` directly (no Python wrapper) — it already
    implements `__iter__` / `__next__` and releases the GIL while blocking.
    """

    def __new__(
        cls,
        ref_path: str,
        vcfs: Mapping[str, str],
        iupac_codes: bool = False,
        missing: str | None = None,
        absent: str | None = None,
        mark_del: str | None = None,
        mark_ins: str | None = None,
        mark_snv: str | None = None,
        mask: str | None = None,
        mask_with: str = "N",
        chain: bool = False,
        regions_overlap: int = 0,
        max_tasks_per_group: int = 0,
        compile_threads: int | None = None,
        log_level: str = "info",
    ):
        return super().__new__(
            cls,
            ref_path=ref_path,
            vcfs=dict(vcfs),
            iupac_codes=iupac_codes,
            missing=missing,
            absent=absent,
            mark_del=mark_del,
            mark_ins=mark_ins,
            mark_snv=mark_snv,
            mask=mask,
            mask_with=mask_with,
            chain=chain,
            regions_overlap=regions_overlap,
            max_tasks_per_group=max_tasks_per_group,
            compile_threads=compile_threads,
            log_level=log_level,
        )

    def __init__(
        self,
        ref_path: str,
        vcfs: Mapping[str, str],
        iupac_codes: bool = False,
        missing: str | None = None,
        absent: str | None = None,
        mark_del: str | None = None,
        mark_ins: str | None = None,
        mark_snv: str | None = None,
        mask: str | None = None,
        mask_with: str = "N",
        chain: bool = False,
        regions_overlap: int = 0,
        max_tasks_per_group: int = 0,
        compile_threads: int | None = None,
        log_level: str = "info",
    ) -> None:
        pass

    # -- convenience: cartesian product -------------------------------------

    def consensus_regions(
        self,
        regions: Sequence[Task],
        samples: Sequence[str] | None = None,
        haplotypes: Sequence[str] | None = None,
        *,
        threads: int = 1,
        prefetch_steps: int | None = None,
        warmup: bool = False,
        ordered: bool = False,
    ) -> _ConsensusIter:
        """Build ``regions × samples × haplotypes`` and return a lazy iterator.

        The cartesian product is materialised eagerly in Python (via
        :func:`build_tasks`) into a flat ``list[Task]``; that list is then handed
        to the Rust engine's lazy iterator verbatim. So "lazy" here means the
        *results* are produced on demand by Rust worker threads — the task list
        itself is fully built up front.

        Pass ``ordered=True`` to get results back in input (region-major) order;
        the default ``ordered=False`` yields in completion order, each result
        carrying its ``idx`` for re-pairing. ``prefetch_steps=None`` lets the
        Rust binding use ``threads`` region groups in flight; pass ``0`` only
        for lowest-memory, one-group-at-a-time iteration.
        """
        tasks = build_tasks(regions, samples, haplotypes)
        return self.consensus_iter(
            tasks,
            prefetch_steps=prefetch_steps,
            warmup=warmup,
            ordered=ordered,
            threads=threads,
        )
