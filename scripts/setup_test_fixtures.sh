#!/usr/bin/env bash
# Set up the external test fixtures for a worktree.
#
# The legacy oracle is a ~400 MB read-only regression fixture that is
# deliberately kept outside git (see docs/test_fixtures.md).  This script
# checks the local fixture copy, restores it from the compressed archive when
# needed, and links it into a worktree as ``.oracle``.
#
# Usage:
#   bash scripts/setup_test_fixtures.sh [WORKTREE]   # default: current dir
#
# Environment:
#   CALIB_FIXTURES_DIR    fixture root, default /home/wqd/cali_data/calib_test_fixtures
#   CALIB_ORACLE_ARCHIVE  archive used when <root>/oracle is missing,
#                         default <root>/oracle_legacy.tar.zst
#   CALIB_ORACLE_URL      optional URL to download the archive from when it is
#                         not present locally
#
# The archive must contain a top-level ``oracle/`` directory and, when present,
# a matching ``<archive>.sha256`` file is verified first.
set -euo pipefail

FIXTURES_DIR="${CALIB_FIXTURES_DIR:-/home/wqd/cali_data/calib_test_fixtures}"
ORACLE_DIR="$FIXTURES_DIR/oracle"
ARCHIVE="${CALIB_ORACLE_ARCHIVE:-$FIXTURES_DIR/oracle_legacy.tar.zst}"
WORKTREE="${1:-$(pwd)}"

verify_archive() {
    local sha="$ARCHIVE.sha256"
    if [ -f "$sha" ]; then
        echo "verifying $(basename "$ARCHIVE")"
        ( cd "$(dirname "$ARCHIVE")" && sha256sum -c "$(basename "$sha")" )
    else
        echo "warning: $sha not found, skipping checksum verification" >&2
    fi
}

if [ ! -d "$ORACLE_DIR" ]; then
    if [ ! -f "$ARCHIVE" ] && [ -n "${CALIB_ORACLE_URL:-}" ]; then
        echo "downloading $ARCHIVE"
        mkdir -p "$FIXTURES_DIR"
        curl -fL --retry 3 -o "$ARCHIVE" "$CALIB_ORACLE_URL"
    fi
    if [ ! -f "$ARCHIVE" ]; then
        echo "error: neither $ORACLE_DIR nor $ARCHIVE exists;" >&2
        echo "       set CALIB_FIXTURES_DIR/CALIB_ORACLE_ARCHIVE/CALIB_ORACLE_URL" >&2
        exit 1
    fi
    verify_archive
    echo "extracting $ARCHIVE into $FIXTURES_DIR"
    mkdir -p "$FIXTURES_DIR"
    tar --zstd -xf "$ARCHIVE" -C "$FIXTURES_DIR"
fi

if [ ! -d "$ORACLE_DIR/09/TB_fit_result" ]; then
    echo "error: $ORACLE_DIR does not look like the oracle fixture" >&2
    exit 1
fi

ln -sfn "$ORACLE_DIR" "$WORKTREE/.oracle"
echo "linked $WORKTREE/.oracle -> $ORACLE_DIR"
echo
echo "raw data must also be mounted per version (just init {ver} {path});"
echo "see docs/test_fixtures.md, or run the oracle test with:"
echo "  CALIB_ORACLE_DIR=$ORACLE_DIR pytest tests/test_pipeline_run.py"