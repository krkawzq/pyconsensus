#!/usr/bin/env python3
"""consensus_enformer.py — drop-in replacement for make_consensus_enformer_new.py.

Reads a genes csv (rows: `gene_id,chr,tss,symbol,strand`) and a samples file
(one sample per line), then produces personal consensus sequences centered on
each gene's TSS using the Rust `pyconsensus` engine instead of spawning
`bcftools consensus` per (gene, sample, haplotype).

Differences from the original script:
  * No intermediate ref fasta slicing — regions are fetched on demand by the
    engine via faidx.
  * No intermediate consensus fasta on disk — sequences are yielded as `bytes`
    by the lazy engine. Here we write them to disk (to mirror the original
    output layout), but a training pipeline can consume them directly.
  * Multi-threaded production inside the engine; the GIL is released while
    workers run.

Usage:
    python consensus_enformer.py <ref_fasta> <genes_csv> <sample_file> -o <out_dir>
    python consensus_enformer.py <ref_fasta> <genes_csv> <sample_file> \\
        --vcf-dir data/variants --threads 16 -o consensus/seqs
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from pyconsensus import ConsensusEngine, Task

# enformer input window (must match the original script / model).
SEQUENCE_LENGTH = 393216
INTERVAL = 114688  # kept for parity with the original script; unused here


def get_vcf_path(vcf_dir: str, chrom: str, pattern: str) -> str:
    """Resolve the VCF for a chromosome.

    `pattern` is a format string with `{chr}` (full, e.g. chr1) and `{n}`
    (number, e.g. 1). Default matches the original GEUVADIS naming.
    """
    name = chrom[len("chr"):] if chrom.startswith("chr") else chrom
    fname = pattern.format(chr=chrom, n=name)
    return os.path.join(vcf_dir, fname)


def load_genes(path: str) -> list[tuple[str, str, int]]:
    """Return [(gene_id, chr, tss)] from a csv with rows gene_id,chr,tss,symbol,strand."""
    genes = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            gene_id, chrom, tss, _symbol, _strand = line.split(",")
            genes.append((gene_id, chrom, int(tss)))
    return genes


def load_samples(path: str) -> list[str]:
    with open(path) as f:
        return [ln.strip() for ln in f if ln.strip()]


def build_regions(genes, vcf_keys):
    """Build `Task` region templates (sample/haplotype left None), TSS-centered.

    1-based inclusive, start clamped to 1 (matches samtools faidx semantics
    and the original script's `if start < 1: start = 1`). The engine expands
    these templates × samples × haplotypes via `consensus_regions`.
    """
    regions = []
    for (gene_id, chrom, tss), vcf_key in zip(genes, vcf_keys):
        start = tss - SEQUENCE_LENGTH // 2
        end = tss + SEQUENCE_LENGTH // 2 - 1
        if start < 1:
            start = 1
        regions.append(Task(chrom, start, end, vcf_key, gene_id))
    return regions


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("ref_fasta", help="reference FASTA (must be faidx-indexed)")
    p.add_argument("genes_csv", help="csv with rows gene_id,chr,tss,symbol,strand")
    p.add_argument("sample_file", help="one sample name per line")
    p.add_argument("--vcf-dir", default="data/variants", help="dir with per-chromosome VCFs")
    p.add_argument(
        "--vcf-pattern",
        default="GEUVADIS.chr{n}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz",
        help="filename pattern with {chr} (full, e.g. chr1) and {n} (number, e.g. 1)",
    )
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("-o", "--out-dir", default="consensus/seqs")
    p.add_argument(
        "--no-write", action="store_true",
        help="don't write fastas (benchmark production only); sequences are still produced",
    )
    args = p.parse_args(argv)

    genes = load_genes(args.genes_csv)
    samples = load_samples(args.sample_file)

    # Map each unique chromosome to its VCF, keyed by chrom name.
    chroms = sorted({chrom for _gid, chrom, _tss in genes})
    vcfs = {}
    for chrom in chroms:
        vpath = get_vcf_path(args.vcf_dir, chrom, args.vcf_pattern)
        if not os.path.exists(vpath):
            sys.exit(f"VCF not found for {chrom}: {vpath}")
        vcfs[chrom] = vpath
    vcf_keys = [chrom for _gid, chrom, _tss in genes]

    regions = build_regions(genes, vcf_keys)

    engine = ConsensusEngine(
        ref_path=args.ref_fasta,
        vcfs=vcfs,
        # match the original script: -I -H {1,2}pIu
        iupac_codes=True,
    )

    # Haplotypes 1pIu / 2pIu. `consensus_regions` expands regions × samples ×
    # haplotypes into a flat task list (Python-side) and feeds it to the Rust
    # lazy iterator verbatim.
    haplotypes = ["1pIu", "2pIu"]
    n_tasks = len(regions) * len(samples) * len(haplotypes)
    print(
        f"Producing {n_tasks} sequences "
        f"({len(regions)} regions × {len(samples)} samples × {len(haplotypes)} haplotypes)",
        file=sys.stderr,
    )

    os.makedirs(args.out_dir, exist_ok=True)

    it = engine.consensus_regions(
        regions,
        samples=samples,
        haplotypes=haplotypes,
        threads=args.threads,
    )

    n_done = 0
    if args.no_write:
        for _idx, _r in it:
            n_done += 1
    else:
        for _idx, r in it:
            sample = r.sample or "_ref"
            sdir = Path(args.out_dir) / sample
            sdir.mkdir(parents=True, exist_ok=True)
            hap = r.haplotype or "all"
            out = sdir / f"{r.gene_id}.{hap}.fa"
            seq = r.seq
            with open(out, "wb") as f:
                f.write(f">{r.gene_id}\n".encode())
                # 60-bp wrap, like bcftools/samtools
                for i in range(0, len(seq), 60):
                    f.write(seq[i : i + 60])
                    f.write(b"\n")
            n_done += 1

    print(f"Done: {n_done} sequences", file=sys.stderr)


if __name__ == "__main__":
    main()
