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
    """Rows: gene_id,chr,tss,symbol,strand -> [(gene_id, chr, tss, symbol, strand)]."""
    genes = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("gene_id"):
                continue
            gid, chrom, tss, symbol, strand = line.split(",")
            genes.append((gid, chrom, int(tss), symbol, strand))
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
def encode_onehot(seq: bytes, include_unknown: bool = True):
    """ACGT (case-insensitive) -> one-hot; anything else -> unknown channel.

    Returns a ``[C, L]`` float tensor:
      * include_unknown=True  -> C=5 (A,C,G,T,N), unknown bases get the N row.
      * include_unknown=False -> C=4, unknown bases are all-zero.

    Lazy import torch so the file is importable without torch installed.
    """
    import torch

    table = b"ACGTacgt"
    n = len(seq)
    channels = 5 if include_unknown else 4
    out = torch.zeros((channels, n), dtype=torch.float32)
    # Vectorised lookup: build a [n] index array, -1 for unknown. bytearray
    # gives a writable buffer for frombuffer without a copy warning.
    arr = torch.frombuffer(bytearray(seq), dtype=torch.uint8)
    idx = torch.full((n,), -1, dtype=torch.long)
    for base in table:  # map both cases of each base to the same channel
        idx[arr == base] = (table.index(base)) % 4
    known = idx >= 0
    if known.any():
        # scatter ones into the right channel for known bases
        ch = idx[known]
        pos = torch.where(known)[0]
        out[ch, pos] = 1.0
    if include_unknown and (~known).any():
        out[4, torch.where(~known)[0]] = 1.0
    return out


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
                yield self._encode(r)

    def _buffer_shuffle(self, it) -> Iterator[dict]:
        import random

        buf: list[dict] = []
        for _idx, r in it:
            buf.append(self._encode(r))
            if len(buf) >= self._shuffle_buffer:
                i = random.randrange(len(buf))
                yield buf.pop(i)
        while buf:
            i = random.randrange(len(buf))
            yield buf.pop(i)

    def _encode(self, r) -> dict:
        return {
            "gene_id": r.gene_id,
            "sample": r.sample,
            "haplotype": r.haplotype,
            "seq": encode_onehot(bytes(r.seq), self._include_unknown),
        }


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
    ):
        self.ref_path = ref_path
        self.vcfs = vcfs
        self.tasks = tasks
        self.threads = threads
        self.shuffle_buffer = shuffle_buffer
        self.include_unknown = include_unknown

    def __iter__(self) -> Iterator[dict]:
        info = _torch_worker_info()
        engine = _worker_state(self.ref_path, self.vcfs)
        # Sharding: each worker takes every info.num_workers-th task, offset by
        # info.id. Strided (not contiguous) sharding keeps gene/sample diversity
        # high within each worker's stream.
        tasks = self.tasks[info.id :: info.num_workers] if info.num_workers > 0 else self.tasks
        it = engine.consensus_iter(tasks, prefetch_steps=16, threads=self.threads)

        # Reuse the single-process encoder + buffer-shuffle on our shard.
        def encode(r) -> dict:
            return {
                "gene_id": r.gene_id,
                "sample": r.sample,
                "haplotype": r.haplotype,
                "seq": encode_onehot(bytes(r.seq), self.include_unknown),
            }

        if self.shuffle_buffer > 1:
            import random

            buf: list[dict] = []
            for _idx, r in it:
                buf.append(encode(r))
                if len(buf) >= self.shuffle_buffer:
                    yield buf.pop(random.randrange(len(buf)))
            while buf:
                yield buf.pop(random.randrange(len(buf)))
        else:
            for _idx, r in it:
                yield encode(r)


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
    p.add_argument("--threads", type=int, default=8, help="Rust worker threads per engine")
    p.add_argument("--num-workers", type=int, default=4, help="PyTorch DataLoader workers (multi mode)")
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--shuffle-buffer", type=int, default=256)
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
            engine, tasks, args.threads, args.shuffle_buffer, include_unknown
        )
        loader = DataLoader(_TorchWrap(ds), batch_size=args.batch_size)
    else:
        # Each worker process builds its own engine from ref_path + vcfs paths;
        # we only resolve the paths here (no engine constructed on the main
        # process for multi mode).
        ds = _MultiWorkerIterableDataset(
            args.ref, _resolve_vcfs(args, genes), tasks, args.threads,
            args.shuffle_buffer, include_unknown,
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
