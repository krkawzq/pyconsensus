#!/usr/bin/env python3
"""enformer_dataloader.py — feed personal consensus sequences into a PyTorch
DataLoader lazily, with zero intermediate fasta on disk.

The pipeline:
  genes (csv) × samples (list) × haplotypes (1pIu, 2pIu)
    -> flat Task list
    -> Rust `consensus_iter` (GIL-released, multi-threaded production)
    -> bytes (393216 bp personal consensus, TSS-centered)
    -> one-hot tensor [4, L] (or [5, L] with an unknown/N channel)
    -> PyTorch DataLoader (multi-process workers + buffer shuffle)

Two composition modes are shown:

  (A) Single-process: one engine, one `consensus_iter`, consumed directly.
      Simplest; the Rust engine already parallelises consensus production
      across `threads` worker threads while Python consumes.

  (B) Multi-process: `DataLoader(num_workers=K)` spawns K worker processes,
      each building its own engine + iterating its shard of the task list.
      Scales to many cores / multi-GPU; ref + VCFs are loaded once per
      worker (the `.cvcf` cache makes subsequent workers' loads cheap).

This file is an *example* — adapt the collate/encoding to your model. It
deliberately avoids writing any fasta to disk; `r.seq` (bytes) goes straight
into a tensor.

Usage:
    # single-process (Rust-threaded)
    python enformer_dataloader.py --ref ref/hg19.fa \\
        --vcf-dir data/Geuvadis/variants/contig_fixed \\
        --genes genes.csv --samples samples.txt --threads 16 --mode single

    # multi-process (PyTorch workers)
    python enformer_dataloader.py --ref ref/hg19.fa \\
        --vcf-dir data/Geuvadis/variants/contig_fixed \\
        --genes genes.csv --samples samples.txt --num-workers 4 --mode multi
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Iterator, Sequence

from pyconsensus import ConsensusEngine, Task, build_tasks

# enformer input window (must match the model + consensus_enformer.py).
SEQUENCE_LENGTH = 393216

# Haplotypes matching the original script (-H {1,2}pIu, -I).
HAPLOTYPES = ("1pIu", "2pIu")

# VCF filename pattern for Geuvadis (per-chromosome, post contig-fix).
VCF_PATTERN = "GEUVADIS.chr{n}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.contig.vcf.gz"


# --------------------------------------------------------------------------- #
# task / region construction
# --------------------------------------------------------------------------- #
def load_genes(path: str) -> list[tuple[str, str, int, str, str]]:
    """Rows: gene_id,chr,tss,symbol,strand -> [(gene_id, chr, tss, symbol, strand)].

    Uses csv parsing so quoted fields (e.g. a symbol containing a comma) and
    ragged rows are handled with a clear, line-numbered error.
    """
    import csv

    genes = []
    with open(path, newline="") as f:
        for lineno, row in enumerate(csv.reader(f), 1):
            if not row or not row[0].strip() or row[0].strip() == "gene_id":
                continue
            if len(row) < 5:
                sys.exit(f"{path}:{lineno}: expected 5 columns, got {len(row)}: {row!r}")
            gid, chrom, tss, symbol, strand = (c.strip() for c in row[:5])
            try:
                tss_int = int(tss)
            except ValueError:
                sys.exit(f"{path}:{lineno}: tss is not an integer: {tss!r}")
            genes.append((gid, chrom, tss_int, symbol, strand))
    return genes


def load_samples(path: str) -> list[str]:
    with open(path) as f:
        return [ln.strip() for ln in f if ln.strip()]


def vcf_path_for(vcf_dir: str, chrom: str, pattern: str) -> str:
    name = chrom[len("chr"):] if chrom.startswith("chr") else chrom
    return os.path.join(vcf_dir, pattern.format(chr=chrom, n=name))


def build_regions(genes: Sequence[tuple[str, str, int, str, str]]) -> list[Task]:
    """TSS-centered region templates (sample/haplotype left None).

    1-based inclusive, start clamped to 1 (samtools faidx semantics; matches
    the original script).
    """
    regions: list[Task] = []
    for gid, chrom, tss, _sym, _strand in genes:
        start = tss - SEQUENCE_LENGTH // 2
        end = tss + SEQUENCE_LENGTH // 2 - 1
        if start < 1:
            start = 1
        # vcf_key == chrom: the engine's vcfs dict is keyed by chrom name.
        regions.append(Task(chrom, start, end, chrom, gid))
    return regions


# --------------------------------------------------------------------------- #
# bytes -> tensor encoding
# --------------------------------------------------------------------------- #
# 256-entry lookup table: byte value -> one-hot channel (0..3 for A/C/G/T),
# or 255 for unknown. Both cases of each base map to the same channel. Built
# once at import as a plain bytearray (no torch dependency); converted to a
# tensor lazily inside encode_onehot so this module imports without torch.
_ENCODER_LUT = bytearray([255]) * 256
for _i, _b in enumerate(b"ACGTacgt"):
    _ENCODER_LUT[_b] = _i % 4


def encode_onehot(seq: bytes, include_unknown: bool = True):
    """ACGT (case-insensitive) -> one-hot; anything else -> unknown channel.

    Returns a ``[C, L]`` float tensor:
      * include_unknown=True  -> C=5 (A,C,G,T,N), unknown bases get the N row.
      * include_unknown=False -> C=4, unknown bases are all-zero.

    One LUT lookup + one scatter per sequence (no per-base Python loop).
    Lazy import torch so the file is importable without torch installed.
    """
    import torch

    n = len(seq)
    channels = 5 if include_unknown else 4
    lut = torch.frombuffer(bytes(_ENCODER_LUT), dtype=torch.uint8)
    arr = torch.frombuffer(bytearray(seq), dtype=torch.uint8)
    ch = lut[arr].to(torch.long)                 # [n], 255 for unknown
    known = ch < 4

    out = torch.zeros((channels, n), dtype=torch.float32)
    if known.any():
        # Scatter 1.0 into [channel, position] for known bases.
        pos = torch.where(known)[0]
        out[ch[known], pos] = 1.0
    if include_unknown and (~known).any():
        out[4, torch.where(~known)[0]] = 1.0
    return out


# Cap on per-process skip warnings; further failures are tallied silently.
_SKIP_WARN_LIMIT = 10
_skip_warn_count = 0


def _encode_result(r, include_unknown: bool) -> dict | None:
    """Turn one ``ConsensusResult`` into a batch dict, or return ``None`` to skip.

    ``r.seq`` is ``bytes | None`` — a failed task (non error-mode) yields
    ``None``. We skip those rather than emit an all-zero tensor that would
    silently poison training, logging the first few with a per-process cap.
    """
    global _skip_warn_count
    if r.seq is None:
        _skip_warn_count += 1
        if _skip_warn_count <= _SKIP_WARN_LIMIT:
            print(
                f"skip {r.gene_id}/{r.sample}/{r.haplotype}: {r.error or 'no sequence'}",
                file=sys.stderr,
            )
        elif _skip_warn_count == _SKIP_WARN_LIMIT + 1:
            print("... further skips suppressed", file=sys.stderr)
        return None
    return {
        "gene_id": r.gene_id,
        "sample": r.sample,
        "haplotype": r.haplotype,
        "seq": encode_onehot(bytes(r.seq), include_unknown),
    }


# --------------------------------------------------------------------------- #
# (A) single-process dataset: one engine, one lazy iterator
# --------------------------------------------------------------------------- #
class ConsensusIterableDataset:
    """Iterates (sample, gene, haplotype) -> one-hot tensors lazily.

    The Rust `consensus_iter` is the producer: it releases the GIL while its
    worker threads produce the next consensus bytes, so Python's consumption
    (encoding) overlaps with Rust's production. No fasta touches disk.

    A small buffer shuffle gives approximate shuffling for SGD without
    materialising the whole (huge) sequence set.
    """

    def __init__(
        self,
        engine: ConsensusEngine,
        tasks: list[Task],
        threads: int,
        shuffle_buffer: int = 256,
        include_unknown: bool = True,
        prefetch_steps: int = 16,
    ):
        self._engine = engine
        self._tasks = tasks
        self._threads = threads
        self._shuffle_buffer = shuffle_buffer
        self._include_unknown = include_unknown
        self._prefetch_steps = prefetch_steps

    def __iter__(self) -> Iterator[dict]:
        it = self._engine.consensus_iter(
            self._tasks,
            prefetch_steps=self._prefetch_steps,
            threads=self._threads,
        )
        if self._shuffle_buffer > 1:
            yield from self._buffer_shuffle(it)
        else:
            for _idx, r in it:
                d = _encode_result(r, self._include_unknown)
                if d is not None:
                    yield d

    def _buffer_shuffle(self, it) -> Iterator[dict]:
        import random

        buf: list[dict] = []
        for _idx, r in it:
            d = _encode_result(r, self._include_unknown)
            if d is not None:
                buf.append(d)
            if len(buf) >= self._shuffle_buffer:
                i = random.randrange(len(buf))
                yield buf.pop(i)
        while buf:
            i = random.randrange(len(buf))
            yield buf.pop(i)


# --------------------------------------------------------------------------- #
# (B) multi-process dataset: each worker owns an engine + a shard
# --------------------------------------------------------------------------- #
def _worker_state(ref_path: str, vcfs: dict[str, str]):
    """Per-worker engine. Built once per worker process; the `.cvcf` cache
    makes the VCF parse cheap after the first worker."""
    if not hasattr(_worker_state, "engine"):
        _worker_state.engine = ConsensusEngine(ref_path=ref_path, vcfs=vcfs, iupac_codes=True)
    return _worker_state.engine


class _MultiWorkerIterableDataset:
    """Shards the task list across DataLoader workers.

    Each worker process builds its own `ConsensusEngine` (ref + VCFs loaded
    once per worker, accelerated by `.cvcf`), takes its slice of tasks, and
    runs `consensus_iter`. `worker_init_fn` seeds per-worker RNG for the
    buffer shuffle.
    """

    def __init__(
        self,
        ref_path: str,
        vcfs: dict[str, str],
        tasks: list[Task],
        threads: int,
        shuffle_buffer: int = 256,
        include_unknown: bool = True,
        prefetch_steps: int = 16,
    ):
        self.ref_path = ref_path
        self.vcfs = vcfs
        self.tasks = tasks
        self.threads = threads
        self.shuffle_buffer = shuffle_buffer
        self.include_unknown = include_unknown
        self.prefetch_steps = prefetch_steps

    def __iter__(self) -> Iterator[dict]:
        info = _torch_worker_info()
        engine = _worker_state(self.ref_path, self.vcfs)
        # Sharding: each worker takes every info.num_workers-th task, offset by
        # info.id. Strided (not contiguous) sharding keeps gene/sample diversity
        # high within each worker's stream.
        tasks = self.tasks[info.id :: info.num_workers] if info.num_workers > 0 else self.tasks
        it = engine.consensus_iter(
            tasks, prefetch_steps=self.prefetch_steps, threads=self.threads
        )

        if self.shuffle_buffer > 1:
            import random

            buf: list[dict] = []
            for _idx, r in it:
                d = _encode_result(r, self.include_unknown)
                if d is not None:
                    buf.append(d)
                if len(buf) >= self.shuffle_buffer:
                    yield buf.pop(random.randrange(len(buf)))
            while buf:
                yield buf.pop(random.randrange(len(buf)))
        else:
            for _idx, r in it:
                d = _encode_result(r, self.include_unknown)
                if d is not None:
                    yield d


class _WorkerInfo:
    __slots__ = ("id", "num_workers")

    def __init__(self, id: int, num_workers: int):
        self.id = id
        self.num_workers = num_workers


def _torch_worker_info() -> _WorkerInfo:
    try:
        from torch.utils.data import get_worker_info
    except ImportError:
        return _WorkerInfo(0, 0)
    info = get_worker_info()
    if info is None:
        return _WorkerInfo(0, 0)
    return _WorkerInfo(info.id, info.num_workers)


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #
def make_engine(ref_path: str, vcf_dir: str, genes, pattern: str) -> ConsensusEngine:
    chroms = sorted({chrom for _gid, chrom, _tss, _s, _st in genes})
    vcfs = {}
    for chrom in chroms:
        p = vcf_path_for(vcf_dir, chrom, pattern)
        if not os.path.exists(p):
            sys.exit(f"VCF not found for {chrom}: {p}")
        vcfs[chrom] = p
    # iupac_codes=True mirrors the original script (-I).
    return ConsensusEngine(ref_path=ref_path, vcfs=vcfs, iupac_codes=True)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ref", required=True, help="reference FASTA (faidx-indexed)")
    p.add_argument("--vcf-dir", required=True, help="dir with per-chrom VCFs")
    p.add_argument("--genes", required=True, help="csv: gene_id,chr,tss,symbol,strand")
    p.add_argument("--samples", required=True, help="one sample name per line")
    p.add_argument("--vcf-pattern", default=VCF_PATTERN)
    p.add_argument("--mode", choices=["single", "multi"], default="single")
    p.add_argument(
        "--threads", type=int, default=8,
        help="Rust worker threads per engine; in multi mode total threads = threads × num_workers",
    )
    p.add_argument("--num-workers", type=int, default=4, help="PyTorch DataLoader workers (multi mode)")
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--shuffle-buffer", type=int, default=256)
    p.add_argument(
        "--prefetch-steps", type=int, default=16,
        help="region groups in flight per engine (memory vs throughput)",
    )
    p.add_argument("--limit", type=int, default=0, help="stop after N batches (0 = all)")
    p.add_argument("--no-unknown-channel", action="store_true", help="4-channel one-hot (no N row)")
    args = p.parse_args(argv)

    genes = load_genes(args.genes)
    samples = load_samples(args.samples)
    regions = build_regions(genes)
    tasks = build_tasks(regions, samples, HAPLOTYPES)
    print(
        f"{len(tasks)} tasks = {len(regions)} genes × {len(samples)} samples × {len(HAPLOTYPES)} haps",
        file=sys.stderr,
    )

    try:
        import torch
        from torch.utils.data import IterableDataset, DataLoader
    except ImportError:
        sys.exit("torch is required for this example: pip install torch")

    # Wrap our dataset in a torch IterableDataset so DataLoader accepts it.
    class _TorchWrap(IterableDataset):
        def __init__(self, ds):
            super().__init__()
            self._ds = ds

        def __iter__(self):
            yield from iter(self._ds)

    include_unknown = not args.no_unknown_channel
    n_done = 0
    if args.mode == "single":
        engine = make_engine(args.ref, args.vcf_dir, genes, args.vcf_pattern)
        ds = ConsensusIterableDataset(
            engine, tasks, args.threads, args.shuffle_buffer, include_unknown,
            prefetch_steps=args.prefetch_steps,
        )
        loader = DataLoader(_TorchWrap(ds), batch_size=args.batch_size)
    else:
        # Each worker process builds its own engine from ref_path + vcfs paths;
        # we only resolve the paths here (no engine constructed on the main
        # process for multi mode).
        if args.threads * args.num_workers > os.cpu_count():
            print(
                f"warning: multi mode spawns {args.threads} threads × {args.num_workers} "
                f"workers = {args.threads * args.num_workers} Rust threads on "
                f"{os.cpu_count()} cores; consider --threads ceil(cores/num_workers)",
                file=sys.stderr,
            )
        ds = _MultiWorkerIterableDataset(
            args.ref, _resolve_vcfs(args, genes), tasks, args.threads,
            args.shuffle_buffer, include_unknown,
            prefetch_steps=args.prefetch_steps,
        )
        loader = DataLoader(_TorchWrap(ds), batch_size=args.batch_size, num_workers=args.num_workers)

    import time

    t0 = time.time()
    for batch in loader:
        # batch is a dict of stacked tensors (seq: [B, C, L]) + string fields.
        seqs = batch["seq"]
        print(
            f"batch: seqs={tuple(seqs.shape)} genes={batch['gene_id']} "
            f"elapsed={time.time() - t0:.1f}s",
            file=sys.stderr,
        )
        n_done += 1
        if args.limit and n_done >= args.limit:
            break

    print(f"done: {n_done} batches in {time.time() - t0:.1f}s", file=sys.stderr)


def _resolve_vcfs(args, genes) -> dict[str, str]:
    chroms = sorted({chrom for _gid, chrom, _tss, _s, _st in genes})
    return {c: vcf_path_for(args.vcf_dir, c, args.vcf_pattern) for c in chroms}


if __name__ == "__main__":
    main()
