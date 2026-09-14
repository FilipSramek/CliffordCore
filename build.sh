#!/usr/bin/env bash
#
# Build and run the CliffordCore test suites.
#
# Every tests/*.cpp is its own translation unit and its own binary -- one per
# algebra, so the algebras never meet in a single TU. Like build_examples.sh
# this keeps going after a failure and reports them all at the end, so one
# run shows every broken suite rather than the first.

set -uo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$ROOT_DIR/build"

mkdir -p "$BUILD_DIR"

# Every header must stand alone. Cheap, and it catches the "works only
# because something else included <string> first" class of bug.
"$ROOT_DIR/check_headers.sh" || exit 1

shopt -s nullglob

total=0
failed=0

for src in "$ROOT_DIR"/tests/*.cpp; do
    name="$(basename "$src" .cpp)"
    total=$((total + 1))

    echo "--- building $name"
    if ! g++ -std=c++17 -I"$ROOT_DIR/include" "$src" -o "$BUILD_DIR/$name"; then
        echo "BUILD FAILED: $name" >&2
        failed=$((failed + 1))
        continue
    fi

    echo "=== $name ==="
    if ! "$BUILD_DIR/$name"; then
        echo "RUN FAILED: $name" >&2
        failed=$((failed + 1))
    fi
done

echo
if [[ $total -eq 0 ]]; then
    echo "No test sources found in $ROOT_DIR/tests/." >&2
    exit 1
fi

echo "$total test suites, $failed failed."
[[ $failed -eq 0 ]] || exit 1
