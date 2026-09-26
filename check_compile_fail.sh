#!/usr/bin/env bash
#
# Every file under tests/compile_fail/ must NOT compile -- and must fail for the
# right reason.
#
# The library rejects non-floating-point component types with a static_assert.
# A static_assert is a hard error, not a substitution failure, so the detection
# idiom the test suites use elsewhere cannot observe it: the only way to test
# that something is rejected is to try compiling it and watch it fail.
#
# A bare "it failed" would prove little, though -- a typo fails too. So each
# file is compiled twice:
#
#   1. with -DCC_T=double, where it must COMPILE. This is the control: it
#      proves the code is otherwise valid.
#   2. as written, with its bad component type, where it must FAIL, and the
#      compiler's output must contain the text on the file's
#      "// expect-error: ..." line.
#
# Run directly, or via ./build.sh which calls it after check_headers.sh.

set -uo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

shopt -s nullglob

total=0
failed=0

for src in "$ROOT_DIR"/tests/compile_fail/*.cpp; do
    name="$(basename "$src")"
    total=$((total + 1))

    expected="$(sed -n 's|^// expect-error: ||p' "$src" | head -n 1 | tr -d '\r')"
    if [[ -z "$expected" ]]; then
        echo "NO EXPECTATION: $name has no '// expect-error: ...' line"
        failed=$((failed + 1))
        continue
    fi

    if ! g++ -std=c++17 -I"$ROOT_DIR/include" -fsyntax-only -DCC_T=double "$src" 2>"$WORK/err.txt"; then
        echo "CONTROL FAILED: $name does not compile even with double"
        sed -n '1,6p' "$WORK/err.txt"
        failed=$((failed + 1))
        continue
    fi

    if g++ -std=c++17 -I"$ROOT_DIR/include" -fsyntax-only "$src" 2>"$WORK/err.txt"; then
        echo "UNEXPECTEDLY COMPILED: $name"
        failed=$((failed + 1))
        continue
    fi

    if ! grep -qF -- "$expected" "$WORK/err.txt"; then
        echo "WRONG ERROR: $name -- expected: $expected"
        sed -n '1,6p' "$WORK/err.txt"
        failed=$((failed + 1))
    fi
done

if [[ $total -eq 0 ]]; then
    echo "No compile-fail tests found under $ROOT_DIR/tests/compile_fail." >&2
    exit 1
fi

echo "$total compile-fail checks, $failed failed."
[[ $failed -eq 0 ]] || exit 1
