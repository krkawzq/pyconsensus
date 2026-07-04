#!/usr/bin/env python3
"""Run a full Geuvadis pyconsensus sink benchmark.

This builds Geuvadis gene x sample x haplotype tasks and calls
``ConsensusEngine.consensus_many_stats`` by default, so consensus sequence bytes
are consumed in Rust and are not returned to Python.

Default paths target:
  /home/wangzhongqi/Code/Project/ExpressionBenchmarkDataset/data/Geuvadis

The current Geuvadis VCF directory in this repo contains a chromosome subset.
By default, genes whose chromosome has no VCF are skipped and reported. Use
``--strict-vcfs`` to fail instead.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

DEFAULT_DATA_DIR = Path(
    "/home/wangzhongqi/Code/Project/ExpressionBenchmarkDataset/data/Geuvadis"
)
DEFAULT_VCF_PATTERN = (
    "GEUVADIS.chr{n}.PH1PH2_465.IMPFRQFILT_BIALLELIC_PH.annotv2.genotypes.vcf.gz"
)
DEFAULT_HAPLOTYPES = ("1pIu", "2pIu")
SEQUENCE_LENGTH = 393_216


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--ref", type=Path, default=None, help="default: DATA_DIR/hg19.fa")
    parser.add_argument(
        "--genes", type=Path, default=None, help="default: DATA_DIR/genes.csv"
    )
    parser.add_argument(
        "--samples", type=Path, default=None, help="default: DATA_DIR/samples.txt"
    )
    parser.add_argument(
        "--vcf-dir",
        type=Path,
        default=None,
        help="default: DATA_DIR/variants",
    )
    parser.add_argument("--vcf-pattern", default=DEFAULT_VCF_PATTERN)
    parser.add_argument("--threads", type=int, default=64)
    parser.add_argument(
        "--compile-threads",
        type=int,
        default=0,
        help="VCF cache/load threads during engine initialization; 0 means auto",
    )
    parser.add_argument(
        "--log-level",
        default="info",
        choices=("off", "error", "warn", "info", "debug"),
        help="Rust engine log level",
    )
    parser.add_argument(
        "--htslib-log-level",
        default="info",
        choices=("off", "error", "warn", "info", "debug", "trace"),
        help="process-global htslib log level",
    )
    parser.add_argument("--haplotypes", nargs="+", default=list(DEFAULT_HAPLOTYPES))
    parser.add_argument("--chrom", action="append", help="restrict to one chromosome; repeatable")
    parser.add_argument("--limit-genes", type=int, default=0, help="0 means no limit")
    parser.add_argument("--limit-samples", type=int, default=0, help="0 means no limit")
    parser.add_argument("--regions-overlap", type=int, default=1)
    parser.add_argument(
        "--max-tasks-per-group",
        type=int,
        default=0,
        help="0 means unlimited grouping, matching engine default",
    )
    parser.add_argument(
        "--strict-vcfs",
        action="store_true",
        help="fail if any gene chromosome has no matching VCF",
    )
    parser.add_argument(
        "--profile",
        action="store_true",
        help="call consensus_many_profile instead of consensus_many_stats",
    )
    parser.add_argument(
        "--compile-stats",
        action="store_true",
        help="print VCF compile stats before the run",
    )
    parser.add_argument(
        "--force-fallback",
        action="store_true",
        help="set PYCONSENSUS_FORCE_FALLBACK_STATE_MACHINE=1 before running",
    )
    parser.add_argument(
        "--disable-same-len-fastpath",
        action="store_true",
        help="set PYCONSENSUS_DISABLE_SAME_LEN_FASTPATH=1 before running",
    )
    parser.add_argument(
        "--disable-edit-fastpath",
        action="store_true",
        help="set PYCONSENSUS_DISABLE_EDIT_FASTPATH=1 before running",
    )
    return parser.parse_args()


def chrom_sort_key(chrom: str) -> tuple[int, str]:
    c = chrom[3:] if chrom.startswith("chr") else chrom
    if c.isdigit():
        return (int(c), "")
    return (10_000, c)


def load_genes(path: Path) -> list[tuple[str, str, int]]:
    genes: list[tuple[str, str, int]] = []
    with path.open(newline="") as handle:
        for row in csv.reader(handle):
            if not row or row[0] == "gene_id":
                continue
            if len(row) < 3:
                raise ValueError(f"bad gene row in {path}: {row!r}")
            genes.append((row[0], row[1], int(row[2])))
    return genes


def load_samples(path: Path) -> list[str]:
    with path.open() as handle:
        return [line.strip() for line in handle if line.strip()]


def vcf_path_for(vcf_dir: Path, pattern: str, chrom: str) -> Path:
    name = chrom[3:] if chrom.startswith("chr") else chrom
    return vcf_dir / pattern.format(chr=chrom, n=name)


def build_vcf_map(
    genes: list[tuple[str, str, int]],
    vcf_dir: Path,
    pattern: str,
    strict: bool,
) -> tuple[dict[str, str], list[tuple[str, str, int]], dict[str, int]]:
    chroms = sorted({chrom for _gid, chrom, _tss in genes}, key=chrom_sort_key)
    vcfs: dict[str, str] = {}
    missing_chroms: list[str] = []
    for chrom in chroms:
        path = vcf_path_for(vcf_dir, pattern, chrom)
        if path.exists():
            vcfs[chrom] = str(path)
        else:
            missing_chroms.append(chrom)

    if missing_chroms and strict:
        missing = ", ".join(missing_chroms)
        raise FileNotFoundError(f"missing VCF(s) for chromosome(s): {missing}")

    missing_set = set(missing_chroms)
    kept_genes = [gene for gene in genes if gene[1] not in missing_set]
    skipped_by_chrom = {
        chrom: sum(1 for _gid, c, _tss in genes if c == chrom) for chrom in missing_chroms
    }
    return vcfs, kept_genes, skipped_by_chrom


def build_tasks(
    genes: list[tuple[str, str, int]],
    samples: list[str],
    haplotypes: list[str],
):
    from pyconsensus import Task

    tasks = []
    half = SEQUENCE_LENGTH // 2
    for gene_id, chrom, tss in genes:
        start = max(1, tss - half)
        end = tss + half - 1
        for sample in samples:
            for hap in haplotypes:
                tasks.append(Task(chrom, start, end, chrom, gene_id, sample, hap))
    return tasks


def print_stderr(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def main() -> int:
    args = parse_args()
    if args.force_fallback:
        os.environ["PYCONSENSUS_FORCE_FALLBACK_STATE_MACHINE"] = "1"
    if args.disable_same_len_fastpath:
        os.environ["PYCONSENSUS_DISABLE_SAME_LEN_FASTPATH"] = "1"
    if args.disable_edit_fastpath:
        os.environ["PYCONSENSUS_DISABLE_EDIT_FASTPATH"] = "1"

    data_dir = args.data_dir
    ref = args.ref or data_dir / "hg19.fa"
    genes_path = args.genes or data_dir / "genes.csv"
    samples_path = args.samples or data_dir / "samples.txt"
    vcf_dir = args.vcf_dir or data_dir / "variants"

    genes = load_genes(genes_path)
    samples = load_samples(samples_path)
    if args.chrom:
        wanted = set(args.chrom)
        genes = [gene for gene in genes if gene[1] in wanted]
    if args.limit_genes > 0:
        genes = genes[: args.limit_genes]
    if args.limit_samples > 0:
        samples = samples[: args.limit_samples]
    vcfs, genes, skipped_by_chrom = build_vcf_map(
        genes, vcf_dir, args.vcf_pattern, args.strict_vcfs
    )

    print_stderr(f"data_dir={data_dir}")
    print_stderr(f"ref={ref}")
    print_stderr(f"genes={genes_path} kept={len(genes)}")
    print_stderr(f"samples={samples_path} n={len(samples)}")
    print_stderr(f"vcf_dir={vcf_dir} n={len(vcfs)}")
    if skipped_by_chrom:
        skipped = sum(skipped_by_chrom.values())
        detail = ", ".join(
            f"{chrom}:{n}" for chrom, n in sorted(skipped_by_chrom.items(), key=lambda x: chrom_sort_key(x[0]))
        )
        print_stderr(f"skipped_genes_without_vcf={skipped} ({detail})")

    tasks = build_tasks(genes, samples, list(args.haplotypes))
    print_stderr(
        "tasks="
        f"{len(tasks)} ({len(genes)} genes x {len(samples)} samples x "
        f"{len(args.haplotypes)} haplotypes)"
    )
    print_stderr(f"threads={args.threads}")
    print_stderr(f"compile_threads={args.compile_threads}")
    print_stderr(f"log_level={args.log_level}")
    print_stderr(f"htslib_log_level={args.htslib_log_level}")

    from pyconsensus import ConsensusEngine, set_htslib_log_level

    set_htslib_log_level(args.htslib_log_level)

    load_started = time.perf_counter()
    engine = ConsensusEngine(
        ref_path=str(ref),
        vcfs=vcfs,
        iupac_codes=True,
        regions_overlap=args.regions_overlap,
        max_tasks_per_group=args.max_tasks_per_group,
        compile_threads=args.compile_threads,
        log_level=args.log_level,
    )
    load_secs = time.perf_counter() - load_started
    print_stderr(f"engine_load_sec={load_secs:.3f}")

    if args.compile_stats:
        for line in engine.compile_stats():
            print(f"compile.{line}")

    run_started = time.perf_counter()
    if args.profile:
        profile_lines = engine.consensus_many_profile(tasks, threads=args.threads)
        run_secs = time.perf_counter() - run_started
        for line in profile_lines:
            print(line)
        result = {
            "mode": "consensus_many_profile",
            "threads": args.threads,
            "tasks": len(tasks),
            "engine_load_sec": load_secs,
            "wall_sec": run_secs,
        }
    else:
        n, total_len, min_len, max_len = engine.consensus_many_stats(
            tasks, threads=args.threads
        )
        run_secs = time.perf_counter() - run_started
        result = {
            "mode": "consensus_many_stats",
            "threads": args.threads,
            "tasks": len(tasks),
            "n": n,
            "total_len": total_len,
            "min_len": min_len,
            "max_len": max_len,
            "engine_load_sec": load_secs,
            "wall_sec": run_secs,
            "seq_per_sec": n / run_secs if run_secs > 0 else None,
            "bases_per_sec": total_len / run_secs if run_secs > 0 else None,
        }

    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
