<div align="center">

# pyconsensus

A Rust rewrite of `bcftools consensus` for the enformer expression-prediction pipeline.

**Personal diploid consensus sequences — online, in-memory, byte-exact.**

[![PyPI](https://img.shields.io/pypi/v/pyconsensus-rs.svg?logo=pypi&logoColor=white)](https://pypi.org/project/pyconsensus-rs/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/Rust-dea584?logo=rust&logoColor=white)](https://www.rust-lang.org/)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Built with PyO3](https://img.shields.io/badge/Built_with-PyO3-orange)](https://pyo3.rs/)

</div>

> Exposed to Python as the native extension package **`pyconsensus`**.

---

## 💡 Why?

The reference pipeline builds a personal diploid consensus sequence for every
*gene × sample × haplotype* combination — on the order of **10M sequences**.
Each sequence is produced by shelling out to `bcftools consensus` as a fresh
subprocess, and the pipeline requires that **all sequences be materialized
upfront, offline, and stored expanded** before downstream code ever runs.

That shape is expensive in three ways that compound across millions of calls:

- **Per-call process startup** — every sequence pays for spawning a process,
  re-opening and re-parsing the *same* VCF, re-reading the *same* reference
  slice, writing a fasta, and building a `.fai`. The consensus algorithm itself
  is cheap; the wrapping I/O and fork overhead is not.
- **Offline pre-materialization** — sequences cannot be produced on demand;
  they must all be generated ahead of time and kept on disk, so the working
  set is the full expanded corpus rather than whatever the consumer currently
  needs.
- **Expanded storage footprint** — every sequence lands on disk as fasta +
  index, so peak storage scales with the total sequence count and length
  rather than with the streaming throughput the consumer can absorb.

The bottleneck is this **spawn-per-sequence + offline + expand-everything**
shape — not the consensus algorithm.

## ⚙️ How

| | Original pipeline | `pyconsensus` |
|---|:---:|:---:|
| Inputs | Reloaded per call | **Loaded once, resident in memory** |
| Consensus | Subprocess per task | Long-lived Rust worker pool |
| Output | Fasta + `.fai` on disk | `bytes` via lazy iterator — **no intermediate files** |
| GIL | — | Released while workers run |

## ✨ Features

- **🔬 Byte-exact parity with `bcftools consensus` 1.23** — a faithful line-by-line port of `forks/bcftools/consensus.c`, covering SNP / MNP / ins / del / complex-indel / `<DEL>` / gVCF `<*>` / `<NON_REF>`, `-H {R,A,I,1,2,1pIu,2pIu,LR,LA,SR,SA}` (+`NpIu`), `-I`, missing/unphased GT, `-a` / `-M` / `--mark-{del,ins,snv}` / `-m --mask-with` / `-c` chain. Pinned by `#[ignore]` tests diffing against a real `bcftools` binary.
- **⚡ Online, no disk** — sequences yielded as `bytes` via a lazy iterator; no intermediate fasta, no temp files.
- **🧠 Resident + multithreaded** — VCFs parsed once into a compiled in-memory store with a binary-searchable region index, allele-operation metadata, compact GT stores, and an optional `.cvcf` cache; `.fai` built once; GIL released while Rust workers run.
- **📈 Built-in profile counters** — `compile_stats()` reports VCF record/allele/GT layout counts, and `consensus_many_profile()` reports lane hit counts, fallback reasons, throughput rates, records seen, and output length totals without returning sequence bytes.
- **📦 Bounded group size** — `max_tasks_per_group=0` keeps legacy unlimited grouping; positive values split large region groups for better load balance and a predictable in-flight task bound.

## Installation

```sh
pip install pyconsensus-rs
```

Pre-built wheels are published on [PyPI](https://pypi.org/project/pyconsensus-rs/)
— no Rust toolchain needed. The wheel is `abi3-py310` and `manylinux_2_34`-compliant,
so a single build covers CPython 3.10+ on x86-64 Linux.

<details>
<summary>Build from source</summary>

```sh
maturin build --release -o dist            # -> dist/pyconsensus_rs-*.whl
pip install dist/pyconsensus_rs-*.whl
```

For local development, `maturin develop --release` installs straight into the
active venv. The first source build downloads htslib plus the compression-lib
tarballs (zlib / bzip2 / xz / libdeflate / zstd) and compiles them with `-fPIC`
so they link cleanly into the cdylib; on restricted networks run `labpon` (or
export `http_proxy` / `https_proxy`) first. `scripts/build_wheels.sh` builds
and verifies the wheel across CPython 3.10–3.15.

</details>

## Usage

```python
from pyconsensus import ConsensusEngine, Task

# The engine preprocesses and holds the inputs (ref + VCFs); the thread count
# is passed per call — the engine itself owns no thread pool.
engine = ConsensusEngine(
    ref_path="ref/hg38.fa",
    vcfs={"chr1": "data/variants/chr1.vcf.gz"},
    iupac_codes=True,          # corresponds to -I
    max_tasks_per_group=8,     # 0 = unlimited; positive values split large region groups
)

# Task coordinates are 1-based inclusive; vcf_key selects an entry of `vcfs`.
# TSS-centered slice (enformer window): start = tss - 393216//2, end = tss + 393216//2 - 1
tss = 2_917_619
start, end = tss - 393216 // 2, tss + 393216 // 2 - 1

tasks = [
    Task("chr1", start, end, "chr1", "ENSG00000263280", "NA12878", "1pIu"),
    Task("chr1", start, end, "chr1", "ENSG00000263280", "NA12878", "2pIu"),
]

# (1) Eager: run a flat task list, get results back in input order.
results = engine.consensus_many(tasks, threads=8)
for r in results:
    print(r.gene_id, r.sample, r.haplotype, len(r.seq or b""))   # r.seq: bytes | None

# (2) Lazy: producer-consumer iterator; GIL released while workers run.
it = engine.consensus_iter(tasks, threads=8, prefetch_steps=16)
for idx, r in it:
    ...   # r.seq feeds the downstream consumer directly

# For very large result streams, batch the Python boundary crossing.
it = engine.consensus_iter(tasks, threads=8, prefetch_steps=16)
while (batch := it.next_batch(256)) is not None:
    for idx, r in batch:
        ...

it = engine.consensus_iter(tasks, threads=8, prefetch_steps=16)
while (batch := it.next_batch_bytes(256)) is not None:
    for idx, seq in batch:
        ...

# (3) Cartesian product (regions × samples × haplotypes) -> lazy iterator.
regions    = [Task("chr1", start, end, "chr1", "ENSG00000263280")]
samples    = ["NA12878", "NA12879"]
haplotypes = ["1pIu", "2pIu"]
for idx, r in engine.consensus_regions(regions, samples, haplotypes, threads=16):
    ...

# (4) Observability: key-value lines suitable for logs.
print("\n".join(engine.compile_stats()))
print("\n".join(engine.consensus_many_profile(tasks, threads=8)))
```

`consensus_iter` / `consensus_regions` yield in **completion order** by default,
each result carrying its input `idx` for re-pairing; pass `ordered=True` to
yield in input order (results are buffered to reassemble). Task expansion is
**region-major** (all haplotypes of a sample for one gene are contiguous, then
the next sample, then the next gene) — matching the original script's layout.

When `max_tasks_per_group > 0`, the number of in-flight tasks is bounded by
`max(prefetch_steps, 1) × max_tasks_per_group`; `0` preserves the original
unlimited group behavior.

## API reference

### `ConsensusEngine(ref_path, vcfs, **opts)`

Constructor options (besides `ref_path` and `vcfs: {vcf_key: path}`):

| Option | Type | Default | bcftools | Notes |
|---|---|---|---|---|
| `iupac_codes` | bool | `False` | `-I` | IUPAC-encode ambiguous alleles |
| `missing` | str\|None | `None` | `-M` | char for missing alleles (single char) |
| `absent` | str\|None | `None` | `-a` | char for absent alleles (single char) |
| `mark_del` | str\|None | `None` | `--mark-del` | char marking deletions (single char) |
| `mark_ins` | str\|None | `None` | `--mark-ins` | `"uc"` / `"lc"` / single char |
| `mark_snv` | str\|None | `None` | `--mark-snv` | `"uc"` / `"lc"` / single char |
| `mask` | str\|None | `None` | `-m` | path to a mask BED |
| `mask_with` | str | `"N"` | `--mask-with` | `"uc"` / `"lc"` / single char |
| `chain` | bool | `False` | `-c` | emit a VCF chain alongside each seq |
| `regions_overlap` | int | `0` | — | record filter: `0`=POS-in-region, `1`=record-span overlap, `2`=variant-span overlap |
| `max_tasks_per_group` | int | `0` | — | `0`=unlimited; positive splits large region groups |
| `compile_threads` | int\|None | `None` | — | rayon pool for VCF compile (defaults to available parallelism capped at the unique-VCF count) |
| `log_level` | str | `"info"` | — | `off` / `error` / `warn` / `info` / `debug` |

### Engine methods

| Method | Returns | Notes |
|---|---|---|
| `consensus_many(tasks, threads=1)` | `list[ConsensusResult]` | input order |
| `consensus_iter(tasks, prefetch_steps=None, warmup=False, ordered=False, threads=1)` | `_ConsensusIter` | lazy; yields `(idx, ConsensusResult)` |
| `consensus_regions(regions, samples=None, haplotypes=None, *, threads=1, prefetch_steps=None, warmup=False, ordered=False)` | `_ConsensusIter` | cartesian product → lazy iterator |
| `consensus_many_stats(tasks, threads=1)` | `(n, total_len, min_len, max_len)` | consumes bytes in Rust; builds no Python list |
| `consensus_iter_stats(tasks, prefetch_steps=None, warmup=False, ordered=False, threads=1)` | `(n, total_len, min_len, max_len)` | iter version of the above |
| `consensus_many_profile(tasks, threads=1)` | `list[str]` | throughput + lane/fallback kv lines |
| `compile_stats()` | `list[str]` | VCF record/allele/GT layout kv lines |
| `log_level` (read/write property) | `str` | runtime engine log level |

### Iterator (`_ConsensusIter`)

| Member | Returns |
|---|---|
| `__iter__` / `__next__` | yields `(idx, ConsensusResult)`; GIL released while blocking |
| `next_batch(batch_size)` | `list[(idx, ConsensusResult)] \| None` |
| `next_batch_bytes(batch_size)` | `list[(idx, bytes)] \| None` — lowest-overhead streaming |

### Module-level functions

```python
from pyconsensus import build_tasks, build_cache, get_htslib_log_level, set_htslib_log_level

# Expand regions × samples × haplotypes into a flat task list (eager).
# A None (or empty) dimension collapses to one task with that field = None.
tasks = build_tasks(regions, samples=["NA12878", "NA12879"], haplotypes=["1pIu", "2pIu"])

# Pre-build .cvcf caches without loading the reference or building the engine.
results = build_cache(
    ["data/variants/chr1.vcf.gz", "data/variants/chr2.vcf.gz"],
    compile_threads=4, force=False,
)
for r in results:
    print(r.path, r.cache_path, r.status, r.records, r.samples, f"{r.cache_mb:.1f}MB")
    # r.status ∈ {"hit", "built", "rebuilt", "forced"}

# htslib log level (off/error/warn/info/debug/trace) — independent of engine log_level.
set_htslib_log_level("error")
print(get_htslib_log_level())
```

### Data classes

| Class | Fields |
|---|---|
| `Task` | `chr`, `start`, `end`, `vcf_key`, `gene_id`, `sample`, `haplotype` (all get/set) |
| `ConsensusResult` | `gene_id`, `sample`, `haplotype`, `seq` (`bytes\|None`), `chain` (`str\|None`), `error` (`str\|None`) — read-only |
| `CacheResult` | `path`, `cache_path`, `status`, `records`, `samples`, `cache_mb`, `elapsed_sec` — read-only |

## License

MIT — see [`LICENSE`](LICENSE).
