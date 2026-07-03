<div align="center">

# pyconsensus

A Rust rewrite of `bcftools consensus` for the enformer expression-prediction pipeline.

**Personal diploid consensus sequences — online, in-memory, byte-exact.**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Rust](https://img.shields.io/badge/Rust-dea584?logo=rust&logoColor=white)](https://www.rust-lang.org/)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Built with PyO3](https://img.shields.io/badge/Built_with-PyO3-orange)](https://pyo3.rs/)

</div>

> Exposed to Python as the native extension package **`pyconsensus`**.

---

## 💡 Why?

The original pipeline (`make_consensus_enformer_new.py`) spawns **one `bcftools consensus` subprocess per *gene × sample × haplotype*** — on the order of **10M invocations**. Each pays for: process start → reopen & reparse the same VCF → re-read the ref slice → write fasta → build `.fai`.

The bottleneck is the **repeated spawning, reparsing, and disk I/O** — not the consensus algorithm itself.

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
- **🧠 Resident + multithreaded** — VCFs parsed once into an in-memory store with a binary-searchable region index; `.fai` built once; GIL released while Rust workers run.

## Build

```sh
maturin build --release -o dist            # -> dist/pyconsensus-*.whl
pip install dist/pyconsensus-*.whl
```

For local development, `maturin develop --release` installs straight into the
active venv. The wheel is `abi3-py310` (abi3), so one build covers CPython 3.10+.

## Usage

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

## License

MIT — see [`LICENSE`](LICENSE).
