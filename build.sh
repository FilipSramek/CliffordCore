#!/usr/bin/env bash
#
# Build and run the CliffordCore test suites.
#
#   ./build.sh              header and compile-fail checks, then every suite
#   ./build.sh pga          the same checks, then only suites whose name contains "pga"
#   ./build.sh --variants   the same, then every selected suite again under each
#                           build variant: -O2, -O2 -DNDEBUG, and -O2 with FMA
#                           contraction (-ffp-contract=fast -mfma, skipped on a
#                           CPU without FMA). Binaries are build/test_<dir>_<variant>.
#
# Every directory under tests/ except common/ and compile_fail/ is one suite --
# one per algebra -- linked into one binary, build/test_<dir>. Each .cpp in it is
# its own translation unit (one per kind of test, plus main.cpp), and they
# compile in parallel; JOBS=<n> overrides the default of one job per CPU. Two
# algebras never meet in one translation unit or one binary. Like
# build_examples.sh this keeps going after a failure and reports them all at the
# end, so one run shows every broken suite rather than the first.
#
# The default build is unoptimised with asserts on. The variants exist because
# that is not what a real-time target ships (docs/safety-audit.md, items A2 and
# B3): -DNDEBUG is the only build that reaches normalize()'s release branch, and
# FMA contraction -- GCC's default for C++ wherever the FPU has it -- changes
# results in the last bits.

set -uo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$ROOT_DIR/build"
JOBS="${JOBS:-$(nproc 2>/dev/null || echo 4)}"
[[ $JOBS =~ ^[1-9][0-9]*$ ]] || JOBS=4

filter=""
variants=0
for arg in "$@"; do
    case "$arg" in
        -h|--help)
            sed -n '2,24p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        --variants) variants=1 ;;
        -*) echo "unknown option: $arg" >&2; exit 2 ;;
        *)  filter="$arg" ;;
    esac
done

mkdir -p "$BUILD_DIR"

# Every header must stand alone. Cheap, and it catches the "works only
# because something else included <string> first" class of bug.
"$ROOT_DIR/check_headers.sh" || exit 1

# Non-floating-point component types must be rejected, with the right message,
# and a discarded [[nodiscard]] result must be diagnosed.
"$ROOT_DIR/check_compile_fail.sh" || exit 1

# Objects are rebuilt from nothing on every run, in a private temporary
# directory: nothing tracks header dependencies, so an old object could never be
# trusted, and a fresh directory means a stale one can never be linked.
WORK="$(mktemp -d)" || exit 1
trap 'rm -rf "$WORK"' EXIT
export CCTEST_INCLUDE="$ROOT_DIR/include"

shopt -s nullglob

total=0
failed=0

# run_suites <variant> [flag...]: build and run every selected suite with the
# given extra flags. An empty variant name is the default build.
run_suites() {
    local variant="$1"
    shift
    export CCTEST_FLAGS="$*"
    local suffix=""
    [[ -n "$variant" ]] && suffix="_$variant"

    for dir in "$ROOT_DIR"/tests/*/; do
        local suite
        suite="$(basename "$dir")"
        case "$suite" in common|compile_fail) continue ;; esac
        [[ -n "$filter" && "$suite" != *"$filter"* ]] && continue
        local name="test_$suite$suffix"
        local objdir="$WORK/$suite$suffix"
        total=$((total + 1))

        local sources=("$dir"*.cpp)
        local objects=()
        mkdir -p "$objdir"
        for src in "${sources[@]}"; do
            objects+=("$objdir/$(basename "$src" .cpp).o")
        done

        echo "--- building $name (${#sources[@]} sources${CCTEST_FLAGS:+, $CCTEST_FLAGS})"
        if [[ ${#sources[@]} -eq 0 ]]; then
            echo "BUILD FAILED: $name has no sources" >&2
            failed=$((failed + 1))
            continue
        fi

        # One g++ per source, at most $JOBS at once. Each job writes its
        # diagnostics to its own log, so parallel output never interleaves, and
        # a failed job leaves no object behind -- the existence check below is
        # the second signal.
        local compiled=1
        for i in "${!sources[@]}"; do
            printf '%s\0%s\0' "${sources[i]}" "${objects[i]}"
        done | xargs -0 -n 2 -P "$JOBS" bash -c \
            'g++ -std=c++17 -Wall -Wextra -Wpedantic $CCTEST_FLAGS -I"$CCTEST_INCLUDE" -c "$1" -o "$2" 2>"$2.log" || { rm -f "$2"; exit 1; }' _ \
            || compiled=0
        for i in "${!sources[@]}"; do
            [[ -s "${objects[i]}.log" ]] && cat "${objects[i]}.log" >&2
            if [[ ! -f "${objects[i]}" ]]; then
                echo "COMPILE FAILED: tests/$suite/$(basename "${sources[i]}")" >&2
                compiled=0
            fi
        done

        # Never leave an old binary behind for a suite that no longer builds.
        rm -f "$BUILD_DIR/$name" "$BUILD_DIR/$name.exe"
        # shellcheck disable=SC2086  # the flags are meant to split
        if [[ $compiled -eq 0 ]] || ! g++ $CCTEST_FLAGS "${objects[@]}" -o "$BUILD_DIR/$name"; then
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
}

run_suites ""

if [[ $variants -eq 1 ]]; then
    run_suites O2 -O2
    run_suites O2_NDEBUG -O2 -DNDEBUG
    if grep -qw fma /proc/cpuinfo 2>/dev/null; then
        run_suites O2_FMA -O2 -ffp-contract=fast -mfma
    else
        echo "--- skipping the FMA variant: this CPU does not report FMA support"
    fi
fi

echo
if [[ $total -eq 0 ]]; then
    if [[ -n "$filter" ]]; then
        echo "No test suites matched '$filter'." >&2
    else
        echo "No test suites found in $ROOT_DIR/tests/." >&2
    fi
    exit 1
fi

echo "$total test suites, $failed failed."
[[ $failed -eq 0 ]] || exit 1
