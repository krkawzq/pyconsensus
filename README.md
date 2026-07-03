# consensus-rs / pyconsensus

A Rust rewrite of `bcftools consensus` that produces personal diploid consensus
sequences **online** for the enformer expression-prediction data pipeline.
Exposed to Python as the native extension package `pyconsensus`.

The original pipeline (`make_consensus_enformer_new.py`) spawns one
`bcftools consensus` subprocess per *gene × sample × haplotype* — on the order
of ten million invocations. Each one pays for: process start → reopen and reparse
the same VCF → re-read the ref slice → write a fasta to disk → build an `.fai`.
The bottleneck is the repeated spawning, reparsing, and disk I/O, not the
consensus algorithm itself.

This tool loads the inputs (reference + VCFs) **once, preprocesses them, and
keeps them resident in memory**; the consensus-apply logic runs in a long-lived
Rust process with a worker thread pool; PyO3 exposes a **lazy iterator** that
produces sequences on demand, yielding `bytes` straight to the training /
inference side — **no intermediate files on disk**.

## Features

- **Byte-exact parity with `bcftools consensus` 1.23** — the apply state machine
  is a faithful line-by-line port of `forks/bcftools/consensus.c`, covering
  SNP / MNP / insertion / deletion / complex-indel / `<DEL>` / gVCF `<*>` /
  `<NON_REF>`, `-H {R,A,I,1,2,1pIu,2pIu,LR,LA,SR,SA}` (and `NpIu`), `-I`,
  missing / unphased GT, `-a` / `-M` / `--mark-{del,ins,snv}` / `-m --mask-with`
  / `-c` chain. Parity is pinned by `#[ignore]` tests that diff against a real
  `bcftools` binary.
- **Online, no disk** — sequences are yielded as `bytes` via a lazy iterator;
  no intermediate fasta, no temp files.
- **Resident + multithreaded** — VCFs are parsed once into an in-memory store
  with a binary-searchable region index; the reference `.fai` is built once; the
  GIL is released while Rust workers run.
- **Fastpath-first execution** — every record is preclassified at load time
  (SNP / same-len / normalized ins / del / symbolic / gVCF …) and each region is
  routed to the cheapest safe lane (in-place patch → edit-script build → full
  state machine), with per-record validation before any output is mutated.
- **Minimal C footprint** — statically links only htslib (+ its bundled
  htscodecs); **bcftools is neither compiled nor linked**, it is an algorithm
  reference only. The sole C glue is a ~30-line struct-field accessor shim.
- **Replaces the script** — ships the TSS-centered slicing logic so it fully
  replaces `make_consensus_enformer_new.py` (see
  [`examples/consensus_enformer.py`](examples/consensus_enformer.py)).

## Build

### Prerequisites

- Rust toolchain (edition 2021), Python ≥ 3.12, [maturin](https://github.com/PyO3/maturin) ≥ 1.0
- C toolchain: `gcc`, `make`, autotools (`autoconf` / `automake` / `libtool`) —
  needed only for htslib's first `autoreconf`
- System decompression libs htslib links against: `zlib`, `bzip2`, `lzma`,
  `zstd` (development headers)

### One-time setup of `forks/` (gitignored)

`forks/` is excluded from version control and must be populated manually:

```
forks/
├── htslib/       # required — built into libhts.a by build.rs
│   └── htscodecs/htscodecs/   # required — htslib's bundled submodule
└── bcftools/     # algorithm reference only — NEVER compiled or linked
```

Place htslib sources under `forks/htslib/` and initialize its bundled htscodecs
submodule:

```sh
cd forks/htslib
git submodule update --init --recursive
```

`build.rs` then auto-runs (only when `libhts.a` is absent):

```sh
autoreconf -i
./configure --disable-libcurl --disable-s3 --disable-gcs --disable-plugins  # CFLAGS=-O2 -fPIC
make -j libhts.a
```

Remote-URL support (curl/s3/gcs) is intentionally disabled — this tool only
reads local files, so the curl dependency is unnecessary. The resulting
`libhts.a` is a minimal static archive linked against the system zlib/bz2/lzma/
zstd. Subsequent builds reuse `libhts.a` unchanged.

### Build the extension

```sh
maturin build --release -o dist            # -> dist/pyconsensus-*.whl
pip install dist/pyconsensus-*.whl
```

For local development, `maturin develop --release` installs straight into the
active venv. The wheel is `abi3-py312`, so one build covers CPython 3.12+.

## Usage

### Python API

```python
from pyconsensus import ConsensusEngine, Task

# The engine preprocesses and holds the inputs (ref + VCFs); the thread count
# is passed per call — the engine itself owns no thread pool.
engine = ConsensusEngine(
    ref_path="ref/hg38.fa",
    vcfs={"chr1": "data/variants/chr1.vcf.gz"},
    iupac_codes=True,          # corresponds to -I
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
    print(r.gene_id, r.sample, r.haplotype, len(r.seq))   # r.seq: bytes

# (2) Lazy: producer-consumer iterator; GIL released while workers run.
for idx, r in engine.consensus_iter(tasks, threads=8, prefetch_steps=16):
    ...   # r.seq feeds the downstream consumer directly

# (3) Cartesian product (regions × samples × haplotypes) -> lazy iterator.
regions    = [Task("chr1", start, end, "chr1", "ENSG00000263280")]
samples    = ["NA12878", "NA12879"]
haplotypes = ["1pIu", "2pIu"]
for idx, r in engine.consensus_regions(regions, samples, haplotypes, threads=16):
    ...
```

`consensus_iter` / `consensus_regions` yield in **completion order** by default,
each result carrying its input `idx` for re-pairing; pass `ordered=True` to
yield in input order (results are buffered to reassemble). Task expansion is
**region-major** (all haplotypes of a sample for one gene are contiguous, then
the next sample, then the next gene) — matching the original script's layout.

### CLI (drop-in for the original script)

[`examples/consensus_enformer.py`](examples/consensus_enformer.py) reads a genes
csv (`gene_id,chr,tss,symbol,strand`) and a samples file (one per line), then
produces TSS-centered consensus sequences via the Rust engine instead of
spawning `bcftools`:

```sh
python examples/consensus_enformer.py \
    ref/hg38.fa genes.csv samples.txt \
    --vcf-dir data/variants --threads 16 -o consensus/seqs
```

Add `--no-write` to produce sequences without writing to disk (for benchmarking
pure production throughput).

## API reference

Public objects (`from pyconsensus import ...`): `ConsensusEngine`, `Task`,
`ConsensusResult`, `build_tasks`, `__version__`. Full signatures live in
[`pyconsensus/_engine.pyi`](pyconsensus/_engine.pyi) (Rust-backed) and
[`pyconsensus/engine.py`](pyconsensus/engine.py) (pure-Python facade).

**`Task(chr, start, end, vcf_key, gene_id, sample=None, haplotype=None)`** — one
production request; **1-based inclusive** coordinates.

**`ConsensusEngine(ref_path, vcfs, **opts)`** — preprocesses and holds the
inputs; the thread count is passed per `consensus_*` call.

| Option | bcftools consensus | Notes |
|---|---|---|
| `iupac_codes` | `-I` | IUPAC codes for heterozygous sites |
| `missing` | `-M` | base used for missing GT (`.`) |
| `absent` | `-a` | base used where no record covers the region |
| `mark_del` / `mark_ins` / `mark_snv` | `--mark-{del,ins,snv}` | `"uc"` / `"lc"` / single char |
| `mask` | `-m` | mask BED file |
| `mask_with` | `--mask-with` | fill for masked regions, default `"N"`; `"uc"`/`"lc"` do not skip variants |
| `chain` | `-c` | also emit the UCSC alignment chain on `ConsensusResult.chain` |
| `regions_overlap` | `--regions-overlap` | `0` = POS in region / `1` = record span / `2` = variant span (default `1`) |
| `Task.haplotype` | `-H` | `R`/`A`/`I`/`1`/`2`/`1pIu`/`2pIu`/`LR`/`LA`/`SR`/`SA` (and `NpIu`) |

**`ConsensusResult`**: `gene_id`, `sample`, `haplotype`, `seq` (`bytes`),
`chain` (`str | None`).

**Engine methods**:

| Method | Returns | Notes |
|---|---|---|
| `consensus_many(tasks, threads=1)` | `list[ConsensusResult]` | parallel; results in input order |
| `consensus_iter(tasks, prefetch_steps=0, warmup=False, ordered=False, threads=1)` | iterator of `(idx, ConsensusResult)` | lazy, GIL-free blocking |
| `consensus_regions(regions, samples=None, haplotypes=None, *, threads=1, prefetch_steps=16, warmup=False, ordered=False)` | iterator of `(idx, ConsensusResult)` | expands the cartesian product first |
| `build_tasks(regions, samples=None, haplotypes=None)` | `list[Task]` | region-major expansion helper |

## How it works

### Layering

```
build.rs              builds forks/htslib into libhts.a (autotools, first build only) + compiles hts_shim.c
src/htslib_ffi.rs     hand-written extern "C" bindings to htslib (faidx / bcf / regidx) — no rust-htslib/noodles
src/hts_shim.c        the only C glue: bcf1_t / bcf_hdr_t field accessors (avoids version-dependent layout)
src/ref_index.rs      faidx wrapper, per-region fetch (Mutex-guarded, Send + Sync)
src/vcf_store.rs      one-shot VCF preprocessing into memory + binary-search region query + .cvcf disk cache
src/compiled.rs       per-record / per-allele preclassification (RecordKind / AlleleOp)
src/planner.rs        region-level fastpath routing (FastPathLane)
src/iupac.rs          IUPAC nucleotide tables, ported from bcftools.h
src/haplotype.rs      -H parsing, sample-mode classification, GT-driven allele selection
src/apply.rs          the consensus apply state machine (port of consensus.c) + fastpath lanes
src/mask.rs           -m mask (regidx-backed) + --mask-with
src/chain.rs          -c UCSC chain output
src/stats.rs          runtime counters (fastpath hits / fallback reasons)
src/engine.rs         ConsensusEngine: task grouping + rayon parallelism + lazy iterator
src/py.rs             PyO3 bindings (feature `python`): Task / ConsensusResult / _ConsensusEngine / _ConsensusIter
pyconsensus/engine.py pure-Python facade: build_tasks + consensus_regions
```

### Coordinate conventions

`faidx_fetch_seq64` uses a **0-based, inclusive, closed** interval; the Python
API and scripts use **1-based inclusive** (matching `samtools faidx chr:start-end`).
The conversion happens at the `RefIndex` boundary, so everything inside the apply
state machine speaks 0-based.

### VCF preprocessing (`vcf_store.rs`)

The VCF is parsed once (eager) into `Vec<VcfRecord>` — pos / rlen / rid /
alleles / per-sample GT / variant type / precompiled metadata. Records are
bucketed by contig and sorted by pos, and a **prefix-max of `ref_end`** is kept
alongside each bucket so a region query can `partition_point`-binary-search the
window **without missing deletions / MNPs whose `pos < start` but whose REF span
reaches into the region**. A `.cvcf` disk cache (keyed on source size + mtime)
avoids re-parsing on repeat runs.

GT decoding goes through `bcf_get_format_values`, so htslib normalizes INT8/16/32
GT width into a flat encoded int32 array; the `bcf_gt_*` bit operations are
replicated in Rust. Biallelic phased diploid records additionally get a bitset
(`BiallelicPhasedGtBits`) so the hot `-H 1/2` path collapses to a single bit
test.

### Fastpath-first execution (`apply.rs` + `planner.rs` + `compiled.rs`)

At load time every record is classified (`RecordKind`: Snp1 / SameLen /
NormInsertion / NormDeletion / SimpleIndel / SymbolicDel / GvcfBlock / Complex)
and every ALT is compiled to an `AlleleOp` (`ref_len` / `alt_len` / `trim_beg` /
`len_diff`). `plan_region` then picks the cheapest safe lane for a region:

1. **Same-len in-place patch** — SNPs / MNPs / same-length replacements applied
   directly to the ref buffer.
2. **Normalized edit-script build** — ins / del / mixed, built cursor-style into
   a fresh buffer with a single `total_delta` length change.
3. **Biallelic phased batch** — for a group of `-H 1/2` tasks over one region,
   a single pass over the records patches every sample's buffer via the GT
   bitset.
4. **Full state machine** — the faithful `apply_variant` port, the correctness
   fallback.

Each fastpath **validates before it mutates** (REF match, bounds, no overlap,
eligible allele op) and falls back to the state machine on any mismatch, so
speed never compromises parity.

### Parallelism (`engine.rs`)

`consensus_many` groups tasks by `(chr, start, end, vcf_key)` so the ref fetch,
VCF query, and region plan are amortized across the sample/haplotype tasks of
one group; identical `(sample, haplotype)` outputs within a group are deduplicated.
`consensus_iter` is a producer-consumer pipeline: a rayon pool (built and torn
down per call) pulls from an unbounded task queue and pushes to a bounded
completion channel (capacity = `prefetch_steps`, giving backpressure); the
iterator blocks on the channel inside `py.allow_threads`, so the GIL is released
while workers produce. `ordered=True` buffers out-of-order arrivals and yields
strictly in input order.

## Tests

```sh
cargo test                     # core-module unit tests (no PyO3, no bcftools needed)
cargo test --features python   # includes the PyO3 binding layer
cargo test -- --ignored        # byte-exact parity vs a real bcftools binary
                               # (needs bcftools / bgzip / tabix on PATH)
```

The `#[ignore]` parity tests (`bcftools_parity`, `bcftools_haplotype_parity`,
`bcftools_m4_parity`) build tiny VCFs + refs on disk, run `bcftools consensus`
with the matching flags, and assert the Rust output is byte-identical.

## Project layout

```
src/            Rust core: apply / haplotype / iupac / mask / chain / vcf_store / compiled / planner / engine / py
build.rs        builds htslib's libhts.a (autotools) on first build, reuses afterwards; compiles hts_shim.c
pyconsensus/    Python package: engine.py (facade) + _engine.abi3.so (Rust ext) + .pyi stubs
examples/       consensus_enformer.py — drop-in replacement for the original bcftools-based script
forks/          (gitignored) htslib / htscodecs / bcftools sources
docs/           (gitignored) design.md / implementation_plan.md / extreme_optimization_plan.md
data/           (gitignored) input datasets
```

## License

MIT — see [`LICENSE`](LICENSE).
