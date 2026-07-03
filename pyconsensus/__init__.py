"""Python package facade for the Rust-backed consensus engine."""

from .engine import ConsensusEngine, build_tasks
from ._engine import ConsensusResult, Task, __version__

__all__ = [
    "ConsensusEngine",
    "ConsensusResult",
    "Task",
    "build_tasks",
    "__version__",
]
