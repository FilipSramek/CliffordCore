#!/usr/bin/env bash
#
# Build and run the CliffordCore examples.
#
#   ./build_examples.sh              build and run everything
#   ./build_examples.sh 03           only examples whose name contains "03"
#   ./build_examples.sh --no-run     build only, do not execute
#
# Binaries go to build/examples/, which is gitignored. Deliberately NOT using
# `set -e`: a broken example should be reported alongside the others, the same
# way the test harness reports every failing check in one run.

set -uo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="$ROOT_DIR/build/examples"

filter=""
run=1
for arg in "$@"; do
    case "$arg" in
        --no-run) run=0 ;;
        -h|--help)
            sed -n '2,12p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        -*) echo "unknown option: $arg" >&2; exit 2 ;;
        *)  filter="$arg" ;;
    esac
done

mkdir -p "$OUT_DIR"
shopt -s nullglob

total=0
failed=0

for src in "$ROOT_DIR"/examples/*.cpp; do
    name="$(basename "$src" .cpp)"
    if [[ -n "$filter" && "$name" != *"$filter"* ]]; then
        continue
    fi
    total=$((total + 1))

    echo "--- building $name"
    if ! g++ -std=c++17 -I"$ROOT_DIR/include" "$src" -o "$OUT_DIR/$name"; then
        echo "BUILD FAILED: $name" >&2
        failed=$((failed + 1))
        continue
    fi

    if [[ $run -eq 1 ]]; then
        echo "=== $name ==="
        if ! "$OUT_DIR/$name"; then
            echo "RUN FAILED: $name" >&2
            failed=$((failed + 1))
        fi
    fi
done

echo
if [[ $total -eq 0 ]]; then
    if [[ -n "$filter" ]]; then
        echo "No examples matched '$filter'." >&2
    else
        echo "No examples found in $ROOT_DIR/examples/." >&2
    fi
    exit 1
fi

echo "$total examples, $failed failed."
[[ $failed -eq 0 ]] || exit 1
