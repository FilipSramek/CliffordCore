#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$ROOT_DIR/build"

mkdir -p "$BUILD_DIR"

g++ -std=c++17 -I"$ROOT_DIR/include" "$ROOT_DIR/tests/test_core.cpp" -o "$BUILD_DIR/test_core"
"$BUILD_DIR/test_core"
