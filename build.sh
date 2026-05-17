#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Configure and build ───────────────────────────────────────────────────────
BUILD_DIR="$SCRIPT_DIR/build"
mkdir -p "$BUILD_DIR"

cmake -S "$SCRIPT_DIR" -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE=Release

cmake --build "$BUILD_DIR" --parallel "$(nproc)"

echo ""
echo "Build complete. Module written to: $SCRIPT_DIR/nd_compress*.so"
