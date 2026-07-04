#!/usr/bin/env bash
#
# build_wheels.sh — build the pyconsensus-rs wheel for CPython 3.10-3.15 and
# verify it installs and imports on every supported version.
#
# The PyO3 binding is compiled with the `abi3-py310` feature (see
# pyproject.toml / Cargo.toml). A single `cp310-abi3` wheel is therefore
# produced that installs on CPython 3.10, 3.11, 3.12, 3.13, 3.14 and 3.15.
# This script builds that wheel once with the lowest interpreter (the abi3
# floor) and then exercises it on each target version, so cross-version
# compatibility is actually verified rather than assumed.
#
# Usage:
#   scripts/build_wheels.sh              # build + verify all versions
#   scripts/build_wheels.sh 3.12         # build + verify a single version
#   scripts/build_wheels.sh 3.11 3.13    # build once (3.11), verify both
#
# Environment:
#   AUDITWHEEL_MODE   "repair" (default) | "warn" | "skip". The compression
#                     libs (z/deflate/bz2/lzma/zstd) are statically linked into
#                     the cdylib (see build.rs), so the .so's NEEDED list is
#                     manylinux-clean and auditwheel repair needs nothing to do.
#   WHEELS_TMPDIR     directory for throw-away verify venvs (default: /tmp).
#
# Prerequisites: uv, maturin, cargo. `uv python install` downloads CPython
# builds and needs outgoing network — on the PJLab dev box run `labpon`
# (or export http_proxy/https_proxy) before invoking this script.

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

VERSIONS=(3.10 3.11 3.12 3.13 3.14 3.15)
if [[ $# -ge 1 ]]; then
  VERSIONS=("$@")
fi

DIST_DIR="$PROJECT_ROOT/dist"
mkdir -p "$DIST_DIR"
: "${AUDITWHEEL_MODE:=repair}"
: "${WHEELS_TMPDIR:=/tmp}"
mkdir -p "$WHEELS_TMPDIR"

# Color-free, grep-friendly progress prefixes. Written to stderr so they never
# pollute values captured via "$(...)" (e.g. interpreter paths).
log()  { printf '>> %s\n' "$*" >&2; }
ver()  { printf '   [%s] %s\n' "$1" "$2" >&2; }

ensure_interpreter() {
  local v="$1"
  local py
  py=$(uv python find "$v" 2>/dev/null || true)
  if [[ -z "$py" ]]; then
    ver "$v" "CPython not found; installing via uv"
    uv python install "$v" >/dev/null
    py=$(uv python find "$v")
  fi
  printf '%s' "$py"
}

# 1. Resolve the build interpreter (lowest requested version = abi3 floor).
BUILD_VER="${VERSIONS[0]}"
log "resolving build interpreter: CPython $BUILD_VER"
BUILD_PY=$(ensure_interpreter "$BUILD_VER")
ver "$BUILD_VER" "build interpreter -> $BUILD_PY"

# 2. Build the release wheel. abi3 produces a single cp310-abi3 artifact
#    regardless of which interpreter builds it; cargo reuses the incremental
#    release cache, so repeated builds are fast.
log "building release wheel (auditwheel=$AUDITWHEEL_MODE)"
maturin build --release \
  --interpreter "$BUILD_PY" \
  --auditwheel "$AUDITWHEEL_MODE" \
  -o "$DIST_DIR"

# 3. Locate the freshly built wheel (only one expected). Match any name so the
#    glob keeps working if the PyPI distribution name ever changes.
mapfile -t WHEELS < <(ls -1 "$DIST_DIR"/*.whl 2>/dev/null || true)
if [[ ${#WHEELS[@]} -eq 0 ]]; then
  log "ERROR: no *.whl found in $DIST_DIR"
  exit 1
fi
WHEEL="${WHEELS[0]}"
log "built: $(basename "$WHEEL")"

# 4. Verify install + import on every target version. Run from /tmp so we
#    load the wheel installed into the venv, not the editable .so that sits
#    next to the source tree in pyconsensus/.
ALL_OK=1
for v in "${VERSIONS[@]}"; do
  log "verifying on CPython $v"
  py=$(ensure_interpreter "$v")
  venv=$(mktemp -d "$WHEELS_TMPDIR/wheel-verify-${v}-XXXXXX")
  uv venv --python "$py" "$venv" >/dev/null
  # UV_LINK_MODE=copy avoids hardlink warnings across filesystems.
  UV_LINK_MODE=copy uv pip install --python "$venv/bin/python" "$WHEEL" >/dev/null 2>&1
  if ( cd /tmp && "$venv/bin/python" -c "
import pyconsensus
from pyconsensus import ConsensusEngine, Task, build_tasks, ConsensusResult, __version__
assert __version__ == '0.1.0'
" 2>/dev/null ); then
    ver "$v" "import OK (__version__=0.1.0)"
  else
    ver "$v" "IMPORT FAILED"
    ALL_OK=0
  fi
  rm -rf "$venv"
done

# 5. Summary.
log "wheels in $DIST_DIR:"
ls -la "$DIST_DIR"/*.whl

if [[ "$ALL_OK" -ne 1 ]]; then
  log "ERROR: one or more versions failed import verification"
  exit 1
fi

log "done. abi3 wheel verified on: ${VERSIONS[*]}"
