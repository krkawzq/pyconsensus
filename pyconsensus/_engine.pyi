from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence

__version__: str

# Runtime accepts:
#   mask_with / mark_ins / mark_snv — "uc", "lc", or any single ASCII char
#       (e.g. "N", the default for mask_with). Typed as plain `str` because a
#       single arbitrary char cannot be expressed as a finite Literal.
#   regions_overlap — 0, 1, or 2 (a u8 on the Rust side). Typed as `int`.
MaskMode: type = str
RegionsOverlap: type = int


class Task:
    """One consensus production request (dataclass-style, Rust-backed)."""

    chr: str
    start: int
    end: int
    vcf_key: str
    gene_id: str
    sample: str | None
    haplotype: str | None

    def __init__(
        self,
        chr: str,
        start: int,
        end: int,
        vcf_key: str,
        gene_id: str,
        sample: str | None = None,
        haplotype: str | None = None,
    ) -> None: ...


class ConsensusResult:
    """One produced consensus sequence (dataclass-style named fields)."""

    gene_id: str
    sample: str | None
    haplotype: str | None
    seq: bytes
    chain: str | None


class _ConsensusIter(Iterator[tuple[int, ConsensusResult]]):
    def __iter__(self) -> _ConsensusIter: ...
    def __next__(self) -> tuple[int, ConsensusResult]: ...
    def next_batch(self, batch_size: int) -> list[tuple[int, ConsensusResult]] | None: ...
    def next_batch_bytes(self, batch_size: int) -> list[tuple[int, bytes]] | None: ...


class _ConsensusEngine:
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
        regions_overlap: int = 1,
        max_tasks_per_group: int = 0,
    ) -> None: ...

    def consensus_many(
        self,
        tasks: Sequence[Task],
        threads: int = 1,
    ) -> list[ConsensusResult]: ...

    def consensus_many_stats(
        self,
        tasks: Sequence[Task],
        threads: int = 1,
    ) -> tuple[int, int, int, int]: ...

    def consensus_many_profile(
        self,
        tasks: Sequence[Task],
        threads: int = 1,
    ) -> list[str]: ...

    def compile_stats(self) -> list[str]: ...

    def consensus_iter(
        self,
        tasks: Sequence[Task],
        prefetch_steps: int | None = None,
        warmup: bool = False,
        ordered: bool = False,
        threads: int = 1,
    ) -> _ConsensusIter: ...

    def consensus_iter_stats(
        self,
        tasks: Sequence[Task],
        prefetch_steps: int | None = None,
        warmup: bool = False,
        ordered: bool = False,
        threads: int = 1,
    ) -> tuple[int, int, int, int]: ...
