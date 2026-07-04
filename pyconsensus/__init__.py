"""Python package facade for the Rust-backed consensus engine."""

from .engine import (
    CacheResult,
    ConsensusEngine,
    build_cache,
    build_tasks,
    get_htslib_log_level,
    set_htslib_log_level,
)
from ._engine import ConsensusResult, Task, __version__

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
