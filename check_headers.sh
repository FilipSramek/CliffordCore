#!/usr/bin/env bash
#
# Every header must compile on its own, with nothing included before it.
#
# This is not busywork: a header that quietly relies on another one having
# already pulled in <cmath> or <string> works fine until someone includes it
# first, and then breaks for a reason that has nothing to do with their code.
# That has happened three times in this repository. The failure is at the
# declaration, not in a template body, so no instantiation is needed to catch
# it -- but only if something actually tries each header alone, which is what
# this script does.
#
# Run directly, or via ./build.sh which calls it first.

set -uo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

total=0
failed=0

while IFS= read -r header; do
    rel="${header#"$ROOT_DIR/include/"}"
    total=$((total + 1))
    printf '#include "%s"\n' "$rel" > "$WORK/probe.cpp"
    if ! g++ -std=c++17 -I"$ROOT_DIR/include" -fsyntax-only "$WORK/probe.cpp" 2>"$WORK/err.txt"; then
        echo "STANDALONE FAILED: $rel"
        sed -n '1,6p' "$WORK/err.txt"
        failed=$((failed + 1))
    fi
done < <(find "$ROOT_DIR/include" -name '*.hpp' | sort)

if [[ $total -eq 0 ]]; then
    echo "No headers found under $ROOT_DIR/include." >&2
    exit 1
fi

echo "$total headers, $failed not standalone."
[[ $failed -eq 0 ]] || exit 1
