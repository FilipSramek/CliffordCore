# Safety-critical and hard real-time readiness audit

**Date:** 2026-10-14 · **Audited:** commit `bd3c8f7`, plus the Phase 1 changes
listed under *What Phase 1 changed* · **Toolchains:** MinGW-W64
g++ 13.2.0 (UCRT), MSVC 14.44 (Build Tools 2022)

This page asks two questions about CliffordCore:

1. **Hard real-time.** Is every call bounded and analysable in time, and free of
   hidden allocation, blocking and shared state?
2. **Safety-critical.** Would the code, and the evidence behind it, survive a
   certification process such as DO-178C (avionics), ISO 26262 (automotive),
   IEC 61508 (industrial) or IEC 62304 (medical)?

No target standard has been chosen yet, so the roadmap uses MISRA C++:2023 as
the coding-standard baseline and notes where the standards differ.

Every code-level finding below was **reproduced**, not just read. The
*Reproductions* section at the end gives the program, the command and the
observed output for each.

## Verdict

| Question | Today | What it would take |
| --- | --- | --- |
| Hard real-time | **Close.** The math core is bounded and analysable by construction. | A handful of code changes (the I/O coupling and the PGA sandwich cost), then measurement on the target |
| Safety-critical | **No.** | Some code changes (error handling above all), but mostly *evidence*: requirements, coverage, a coding standard, a qualified toolchain |

The foundations are unusually good for real-time use. The gap to safety use is
mostly process rather than code, and the code part is concentrated in one
place: what the library does with degenerate input.

## What is already right

- The math path has **no heap allocation, no exceptions, no RTTI, no virtual
  functions, no loops, no recursion, no `goto` and no global or static mutable
  state**. Every function runs a fixed number of operations with a few
  data-dependent branches, so its worst-case execution time (WCET) is bounded.
  Every function is reentrant, and it is thread-safe on distinct objects.
- Default constructors zero every member, so nothing is ever left uninitialized.
- `static_assert(is_floating_point)` rejects integer and `bool` component types,
  and `tests/compile_fail/` pins the rejection.
- The tests use **independent oracles**: the blade-product rule for the Cayley
  tables, and `std::complex` for Cl2. Property tests run in float, double and
  long double. That is the right kind of verification evidence.
- The numerics are careful: `atan2` rather than `acos`, and a `rotor_between`
  that does not cancel catastrophically near opposite directions.
- It is header-only with no third-party dependency, which keeps the
  certification surface small: about 17k header lines and 1,100 function
  definitions.
- It is C++17, which is the language MISRA C++:2023 targets.
- **NaN is contagious everywhere it should be.** The new sweep puts a NaN in
  every component of every operand of the public entry points, in all three
  floating-point types. Only the PGA functions that ignore `e0` by design
  produce a finite result (see A9). Nothing turns a NaN into a plausible value.
- **The suites pass under every build tried:**
  - g++ at `-O0`, `-O2`, `-O2 -DNDEBUG` and `-O2 -ffp-contract=fast -mfma`
  - MSVC `/W4 /permissive-` in debug and in `/O2 /DNDEBUG`

## Findings

Severity: **B** = blocker for any safety use · **M** = major · **m** = minor.
Status: **Open**, **Fixed** (in Phase 1), or **Partial** (Phase 1 reduced it; the
rest is listed).

### A. Error handling and numerics

| # | Sev | Status | Finding |
| --- | --- | --- | --- |
| A1 | B | Open | **A degenerate input returns a plausible-looking value**, so the caller cannot tell "invalid" from a real answer. `rotor_from_axis_angle(zero axis)` and `rotor_between(zero, v)` give the identity rotor, and PGA `motor_between(ideal plane, p)` the identity motor. PGA `distance(point, plane)` gives `0` for an ideal point or plane, which reads as "on the plane". `position(ideal point)` gives the origin. Cl2 `angle(zero)` and `angle_between(zero, v)` give `0`. `log(zero rotor)` gives the zero bivector, which is also `log` of the identity. A release-build `normalize(0)` returns its input, which is not unit length. All of these are pinned by `special_values.cpp`, so changing them is deliberate. |
| A2 | B | Partial | **Debug and release builds behave differently.** `normalize` asserts (aborts) when asserts are on, and returns its input unchanged under `NDEBUG`. Until Phase 1 the release branch never ran in any test. `./build.sh --variants` now runs it, via `test_normalize_zero_release`, but the divergence itself remains: what is tested in debug is still not what ships. |
| A3 | B | Open | **The algebras handle division by zero inconsistently.** Cl3 and Cl2 single-grade `inverse(0)` gives Inf or NaN (as do `project`, `reflect` and `reject`, which call it). PGA's give zero, and `inverse(Multivector)` gives zero in every algebra. MSVC at `/O2` flags the division itself (warning C4723) when it can see a constant zero. |
| A4 | B | Partial | **There is no overflow-safe norm.** Every `norm` is a naive `sqrt(x² + y² + z²)`. Outside the *Input domain* below, results are silently wrong. For example, `rotor_from_axis_angle(Vector<float>(1e20, 0, 0), π/2)` returns the identity rotor, and PGA `inverse(Bivector)` drops its moment part once a float line's direction reaches 2^32. Phase 1 defined the domain and tests both of its edges. The norms themselves are unchanged. |
| A5 | B | **Fixed** | **The PGA projection bug.** `project(Line, Plane)` returned the line reversed, and the test hid it with `std::abs`. It now uses the right contraction, `(L ⌊ p) p`, and the test checks orientation against the line through the two projected points. |
| A6 | M | Partial | **Preconditions on rotors and motors are unchecked.** `sandwich`/`rotate`/`transform`/`translate` assume a unit versor: a rotor that has drifted to scalar 1.1 scales `e1` to `1.21 e1`. `log(Motor)` keeps the rotation under drift but not the translation, which 1% drift moves by about 1%. All 38 overloads and `log(Motor)` now document this with `@pre`, but nothing checks it, and nothing in the type system separates a unit rotor from any other. |
| A7 | m | Open | **`approx_equal` uses an absolute tolerance only** (100 ε). It is too strict for components around 1e4 and meaningless around 1e-14. That is harmless in tests and wrong for run-time decisions. |
| A8 | M | Open | **In debug, `normalize` aborts on a non-zero input**, with the message "the vector is zero". A vector of 1e-200 has a squared norm that underflows to zero. In release, the same call returns the input unchanged, which is not unit length. |
| A9 | M | Open | **PGA `norm`, `angle` and the inner product `\|` ignore every `e0` component by design**, so a NaN in an ideal component is invisible to them. `norm(Vector(NaN, 1, 0, 0))` is exactly 1. A sanity check of the form `isfinite(norm(x))` therefore passes on a poisoned plane, and `is_ideal` of a point with a NaN weight is `false`. This is pinned by `test_norm_ignores_ideal_nan`. |

### B. Real-time behaviour

| # | Sev | Status | Finding |
| --- | --- | --- | --- |
| B1 | M | Open | **Every type header includes `<string>`, and through `detail/format.hpp` `<sstream>`, because `to_string()` is a class member.** That brings iostreams, locales, the heap and exceptions into every translation unit's include set. It blocks builds on bare-metal or restricted C++ libraries, and conflicts with MISRA C++:2023 Rule 21.6.1 and JSF AV Rule 206 on dynamic memory. Uncalled templates generate no code, but the dependency is still there. |
| B2 | M | Open | **The PGA sandwich runs the full 16×16 table twice.** Measured on the host at `-O2`, `sandwich(point, motor)` takes **273 ns**, against **7.7 ns** for Cl3's closed-form `rotate(vector, rotor)`, about 35×. `sandwich(point, rotor)` takes 302 ns, because it widens to a motor first. IEEE semantics forbid deleting the `x * 0.0` terms, so they all execute at any optimization level. It is the hottest operation in a rigid-body loop. |
| B3 | M | Partial | **Optimized and contracted builds give different results.** The default build is `-O0`. GCC's C++ default is `-ffp-contract=fast`, which fuses multiply-adds wherever the FPU has FMA: AArch64, Cortex-M4F/M7 and x86 with `-mfma`. With contraction on, `a ^ a` is about 1e-17 rather than exactly 0, and a transformed point's weight is 0.99999999999999989 rather than 1. Host and target can therefore disagree in the last bits, and an exact-zero guard can miss. The suites pass under contraction, and `--variants` now runs that build. Testing on an actual target remains open. |
| B4 | M | Open | **`sin`, `cos`, `atan2` and `sqrt` take their timing from the target's math library.** They are used by exp/log/slerp, rotor and motor construction, `distance` and `angle`. A very large angle, such as 1e6 rad, can hit slow argument reduction. Bound the input angle range, and measure on the target. |
| B5 | m | Open | **Subnormal numbers can cost 100+ cycles per operation** on x86 and many other FPUs, and squared norms of small components produce them. Stay inside the input domain, or run the target with flush-to-zero. |

### C. Coding-standard conformance

The rule numbers are given for orientation. Confirm them against the text of the
standard you adopt.

| # | Sev | Status | Finding | Rule |
| --- | --- | --- | --- | --- |
| C1 | M | Open | **Single-argument constructors are implicit:** `Scalar(T)`, `Trivector(const T&)`, `Quadvector(const T&)`, Cl2 `Bivector(const T&)`, and the widening `Multivector(Rotor)`, `Motor(Rotor)` and `Motor(Translator)`. The first group is also the cause of the known `Scalar * T` ambiguity. | MISRA C++:2023 15.1.3; AUTOSAR A12-1-4 |
| C2 | M | **Fixed** | **Discarded results compiled silently.** `normalize(r);`, written by someone who expected it to normalize in place, compiled with no diagnostic in g++ (`-Wall -Wextra -Wpedantic`) or MSVC (`/W4`). All 962 functions that return a value are now `[[nodiscard]]`, and three compile-fail checks pin it. | MISRA C++:2023 0.1.2 |
| C3 | M | Open | **Floating-point values are compared with exact `==`**, in guards (`n == T(0)`) and in the public `operator==`. The guards miss near-zero values; see A8. | MISRA C++:2008 6-2-2; AUTOSAR M6-2-2 |
| C4 | M | Open | **Functions marked `constexpr` call `std::sqrt`, `sin`, `cos`, `atan2` and `abs`, none of which is `constexpr` in C++17 or C++20.** g++ accepts `constexpr auto n = norm(v);` because it treats these as builtins. MSVC rejects it (error C2131, at `sqrt`). The code is portable at run time, but the compile-time claim is GCC-only, and formally ill-formed with no diagnostic required. | ISO C++17 [dcl.constexpr] |
| C5 | m | Open | **`#include` appears after code** (the two-phase declare-then-define pattern). It is deliberate and load-bearing, so record it as a documented deviation rather than restructuring. | MISRA C++:2008 16-0-1; MISRA C++:2023 19.0.3 |
| C6 | m | Open | **`^` and `\|` bind more loosely than `+`, `-`, `*` and `==`**, so `a ^ b + c` means `a ^ (b + c)`. The types catch most misuse, but not when the result goes into `auto`. Named forms exist (`wedge_product`, `left_contraction`); a project guideline should require them, or parentheses. | MISRA C++:2023 8.0.1 |
| C7 | m | **Fixed** | **No function was `noexcept`.** All 1,111 functions and constructors now are, except the printing ones (`to_string`, `operator<<`, `format_component`), which allocate. | AUTOSAR A15-4-x |
| C8 | m | Open | **`format_component` follows the global locale**, so a de-DE locale prints `1,5`. | — |

### D. Process and evidence

| # | Sev | Status | Gap |
| --- | --- | --- | --- |
| D1 | B | Open | There are no requirements and no traceability from tests to requirements. The existing tests are good raw material for both. |
| D2 | B | Open | No structural coverage is measured. See *What the standards add* for the level each one requires. For templates, coverage counts per instantiation actually used. |
| D3 | B | Partial | There is no static analysis: no MISRA checker, clang-tidy or cppcheck. A one-off UBSan run in trap mode (MinGW ships no sanitizer runtime) found no undefined behaviour, and also showed that **no test reached a division by zero** before Phase 1. The new zero-input tests now do so deliberately, which is why trap-mode `float-divide-by-zero` cannot become a permanent build variant. |
| D4 | M | Open | Two toolchains have been tested, both on a Windows x86-64 host. There is no CI, no ARM cross-build (`arm-none-eabi-g++ -fno-exceptions -fno-rtti`), and no trial with a qualified compiler such as IAR, Green Hills, TASKING or Arm FuSa. |
| D5 | M | Open | Configuration management: there are no release tags, and the version disagrees between Cl2 (0.3.0) and the rest (0.2.0). Certification needs fixed, identifiable baselines. |
| D6 | M | Partial | There is no safety manual. In ISO 26262 terms the library would be a Safety Element out of Context (SEooC), with assumed requirements and assumptions of use. *Assumptions of use* below is the seed of one. |
| D7 | m | Open | There is one developer. The higher integrity levels require independent verification. |

## Input domain

**Every component of every input, intermediate and output should have magnitude
in [2^-30, 2^30] in float, and in [2^-250, 2^250] in double and long double.**
Zero is always allowed. Inside that range no expression in the library overflows
or underflows. Outside it, some results are silently wrong rather than flagged.

The bound comes from the tightest expression in the library, `a * a` in PGA
`inverse(Bivector)`, where `a` is itself a squared norm. In float, that line's
moment part survives at a direction of 2^30 and is lost at 2^32. Most functions
tolerate far more. The range is a sufficient condition, not a per-function
limit. In physical terms, float in metres covers about a nanometre to a million
kilometres, but **products grow magnitudes**. Joining two points about 2^30
apart already produces a line at the edge. Check the domain for chained
expressions, not just raw inputs.

Long double uses the double range: it is 80-bit on x86 but the same as double
on ARM and with MSVC, so the double range is the portable one.
`special_values.cpp` tests accuracy at both edges, for all three types and all
three algebras.

## Assumptions of use

These are the conditions under which the tested behaviour holds. A safety case
built on this library would carry them as requirements on the integrator.

- **Floating point.** IEEE 754 binary32 or binary64, round-to-nearest, and **no
  `-ffast-math`** (and no `-ffinite-math-only`). Those flags let the compiler
  assume NaN never occurs, which removes exactly the contagion guarantee above.
- **FMA contraction.** It is allowed, and the suites pass under it, but results
  then differ from an uncontracted build in the last bits. If host and target
  must agree bit for bit (replay testing, or comparing redundant channels),
  compile with `-ffp-contract=off`. GCC ignores `#pragma STDC FP_CONTRACT` in
  C++, so it has to be a flag.
- **Inputs inside the domain** above. Check the domain at the boundary of the
  system, where data enters.
- **Unit versors.** Every rotor, motor and translator passed to `sandwich`,
  `rotate`, `transform`, `translate`, `slerp` or `log(Motor)` must be unit
  (scalar 1 for a translator). Renormalize after composition. Rotors drift
  faster in float than in double.
- **No reliance on degenerate-input results.** Do not pass zero or ideal
  elements where a direction, a weight or an inverse is needed, and do not
  treat what comes back as an error signal. The table below records what each
  returns today.
- **No reliance on asserts.** Asserts are a debugging aid only: debug aborts,
  release returns the input unchanged.
- **NaN checks on PGA ideal parts are separate.** `norm` does not see them
  (A9). Check `std::isfinite` on the components themselves.
- **Printing stays out of real-time paths.** `to_string`, `operator<<` and
  `format_component` allocate and can throw. Everything else is `noexcept`.
- **Compile-time evaluation** of anything that calls the math library works
  only on GCC (C4). Do not rely on it.
- **Subnormals.** Either stay inside the domain or run with flush-to-zero, if
  timing matters more than gradual underflow.

### What degenerate inputs return today

These cases are pinned by `tests/<algebra>/special_values.cpp` and
`test_normalize_zero_release`:

| Call | Cl2 | Cl3 | PGA |
| --- | --- | --- | --- |
| `inverse(Scalar 0)` | +Inf | +Inf | +Inf |
| `inverse` of a zero or ideal blade, rotor or translator | NaN | NaN | zero |
| `inverse(Multivector 0)` | zero | zero | zero |
| `normalize(0)` or `normalize(ideal)`, debug | abort | abort | abort |
| `normalize(0)` or `normalize(ideal)`, release | input unchanged | input unchanged | input unchanged |
| `rotor_between(0, v)`, `rotor_between(v, 0)` | identity | identity | — |
| `rotor_from_axis_angle(zero axis, a)` | — | identity | identity |
| `motor_between(ideal plane, p)` | — | — | identity |
| `log(zero rotor)` or `log(zero motor)` | zero bivector | zero bivector | zero line |
| `angle(0)`, `angle_between(0, v)`, `rotor_angle(0)` | 0 | — | — |
| `distance` involving an ideal point, plane or line | — | — | 0 |
| `position(ideal point)` | — | — | origin |
| `project`, `reflect`, `reject` onto a zero direction | NaN | NaN | zero |

## Real-time characteristics

- **Bounded by construction.** There are no loops and no recursion, so each
  call's op count is fixed. The data-dependent branches are a handful of
  zero tests and sign tests (`rotor_between`, `slerp`), and each takes one of
  two fixed paths.
- **No allocation, no exceptions, no blocking** in anything except the printing
  functions.
- **Stack use is small and fixed.** The largest value type is PGA's
  `Multivector<double>`, at 128 bytes. A sandwich holds a few such temporaries.
  Measure with `-fstack-usage` on the target for the budget.
- **Measured costs** on the host (x86-64, g++ 13.2 `-O2`, a dependency chain of
  2,000,000 calls, so latency rather than throughput):

  | Operation | ns per call |
  | --- | --- |
  | Cl3 `rotate(Vector, Rotor)` (closed form) | 7.7 |
  | PGA `sandwich(Trivector, Motor)` (two 16×16 products) | 273 |
  | PGA `sandwich(Trivector, Rotor)` (widened to a motor) | 302 |

  These numbers are relative evidence only. WCET has to be measured, or
  analysed statically, on the target hardware with its own math library.
- **Transcendental calls**, which are the only externally timed code, appear
  in:
  - `exp`, `log` and `slerp`
  - rotor and motor construction (`rotor_from_axis_angle`, `rotor_from_angle`,
    `motor_from_line_angle`, `screw`)
  - `rotor_angle`, `angle`, `angle_between` and `distance`
  - every `norm`, through `sqrt`

## What Phase 1 changed

Phase 1 is the set of non-breaking changes made alongside this audit:

| Change | Where |
| --- | --- |
| A5 fixed: `project(Line, Plane)` keeps the line's orientation; the masking test was replaced by signed checks | `include/cliffordcore/pga/operations/geometry.hpp`, `tests/pga/geometry.cpp` |
| C2: `[[nodiscard]]` on all 962 functions that return a value; three compile-fail checks (`cl2_discarded_sum`, `cl3_discarded_normalize`, `pga_discarded_inverse`) | every header; `tests/compile_fail/` |
| C7: `noexcept` on all 1,111 functions and constructors except the printing ones | every header |
| A2: `test_normalize_zero_release` in each suite, compiled only under `NDEBUG` | `tests/*/unary_ops.cpp` |
| A1/A3/A4/A9: new `special_values` category: the NaN-contagion sweep (float, double, long double), zero and ideal inputs pinned, accuracy at both domain edges, and the PGA `e0` blind spot | `tests/*/special_values.cpp` |
| B3: `./build.sh --variants` reruns every suite under `-O2`, `-O2 -DNDEBUG` and `-O2 -ffp-contract=fast -mfma` | `build.sh` |
| A6: a unit-versor `@pre` on all 38 `sandwich`/`rotate`/`transform`/`translate` overloads and on `log(Motor)` | `*/operations/sandwich.hpp`, `pga/operations/log.hpp` |

Check counts after Phase 1:

| Build | Cl2 | Cl3 | PGA |
| --- | --- | --- | --- |
| default, `-O2`, `-O2` with FMA | 4519 | 3357 | 4930 |
| `-O2 -DNDEBUG` | 4525 | 3364 | 4939 |

## Roadmap

### Phase 2: API decisions (weeks; breaking)

1. **Choose one degenerate-input policy** and make both builds behave the same
   (A1–A3, A8). The cheapest form that fits hard real-time is a checked layer
   beside the unchecked one: `try_normalize`, `try_inverse` and so on, each
   returning `{value, status}`, with no heap and no exceptions. The existing
   functions stay the fast path for inputs already validated, with their
   preconditions documented.
2. **Add `is_finite(x)` for every type**, which closes A9's blind spot without
   changing `norm`.
3. **Move `to_string` out of the classes**, into free functions in the existing
   `operations/stream.hpp`. The math headers would then include only `<cmath>`,
   `<type_traits>`, `<limits>` and `<array>` (B1).
4. **Make the constructors `explicit`** (C1). This also fixes `Scalar * T`.
5. **Add closed-form PGA sandwiches** for points, planes and lines, checked
   against the general product the same way `Rotor * Rotor` already is (B2).
6. **Handle overflow:** use overflow-safe norms (scaled, or `std::hypot`, which
   is slower), or enforce the domain at entry points in the checked layer (A4).
7. **Settle `constexpr`** on functions that call the math library (C4). Either
   document it as GCC-only, or remove it from those functions.
8. **Use relative-plus-absolute tolerance** in `approx_equal` (A7), and
   **format numbers locale-independently** in `format_component` (C8).

### Phase 3: certification evidence (months; needed regardless of code quality)

1. **Choose the target standard and integrity level.** Everything below scales
   with it.
2. **Requirements and traceability:** write a requirements specification and a
   trace matrix from requirements to tests (D1).
3. **Coverage:** measure statement and branch coverage with gcov/lcov, and
   MC/DC with a dedicated tool where the level requires it (D2).
4. **Static analysis:** run a MISRA C++:2023 checker, with recorded deviations
   for C5 and for operator overloading, plus clang-tidy or cppcheck in CI (D3).
5. **CI** across several compilers, including an ARM cross-build with
   `-fno-exceptions -fno-rtti`, ideally plus a qualified compiler (D4).
6. **Release baselines:** tagged releases and consistent versions (D5).
7. **Safety manual:** the assumptions of use above, completed (D6).
8. **WCET:** measurement on target hardware, or a static WCET tool (B4).
9. **Independent review** of the requirements, the tests and the code (D7).

### What the standards add

| Topic | DO-178C | ISO 26262 | IEC 61508 |
| --- | --- | --- | --- |
| Integrity levels | DAL A (highest) to E | ASIL D (highest) to A | SIL 4 (highest) to 1 |
| Structural coverage | statement (DAL C), plus decision (DAL B), plus MC/DC (DAL A) | statement and branch rising to MC/DC, highly recommended at ASIL D | rising with SIL, to MC/DC at the top |
| Library route | a previously developed component, verified like the rest of the software | SEooC, or qualification of a software component (Part 8) | assessed as part of the software, with a language subset |
| Coding standard | a project standard is required; MISRA C++ is the usual choice | a language subset is required; MISRA C++:2023 absorbed AUTOSAR C++14 | a language subset is required; MISRA C++:2023 fits |
| Tools | qualification (DO-330) for any tool whose output is not verified | tool confidence level (Part 8) | tool classes T1–T3 |

## Reproductions

Each finding was reproduced with a throwaway program compiled against
`include/`. The programs are not part of the repository; the snippets below
are the parts that matter.

**A1, A3, A4, A6, A8 (Cl3).** Built with
`g++ -std=c++17 -Iinclude probe.cpp`, and again with `-O2 -DNDEBUG`.

```text
rotor_from_axis_angle(zero axis, 1.0)                 s=1 xy=0 xz=0 yz=0
rotor_between(zero, e1)                               s=1 xy=0 xz=0 yz=0
inverse(Vector 0)                                     (nan, nan, nan)
inverse(Scalar 0)                                     inf
inverse(Multivector 0)                                scalar=0
rotor_from_axis_angle(Vector<float>(1e20,0,0), pi/2)  s=1 xy=0 xz=0 yz=0        <- identity
rotor_from_axis_angle(Vector<float>(1,0,0),   pi/2)   s=0.707107 ... yz=-0.707107
normalize(Vector<float>(1e20,0,0))                    (0, 0, 0)
rotate(e2, rotor_between(Vector<float>(3e19,0,0), e2))  (0, 1, 0)                <- should be (-1, 0, 0)
norm(Vector(1e-200,0,0))                              0
rotate(e1, rotor with scalar 1.1)                     (1.21, 0, 0)
release:  normalize(Vector(1e-200,0,0))               (1e-200, 0, 0)            <- not unit
debug:    Assertion failed: n != T(0) && "normalize: the vector is zero and has no unit direction"
```

**A5 (PGA, before the fix).** The line from (0, -2, 0) to (1, -2, 0) was
projected onto z = 1:

```text
direction(line)                                (1, 0, 0)
direction(project(line, z=1))                  (-1, -0, 0)      <- reversed
direction(join(project(p1), project(p2)))      (1, 0, 0)
direction(grade2((L |_ p) p) / |p|^2)          (1, -0, 0)       <- the fix
tilted, onto y = 3: current (-0.707,0,-0.707), via points (0.707,0,0.707), fixed (0.707,0,0.707)
```

**A1, A4, A9 (PGA).**

```text
distance(ideal_point(1,0,0), plane z=1)        0
distance(point(0,0,5), ideal_plane)            0
position(ideal_point(1,2,3))                   (0, 0, 0)
inverse(Vector 0)                              (0, 0, 0, 0)
inverse(Bivector(1, 0,0, 2^k, 0,0)) in float:  k=30  e01=-8.67e-19 (correct)
                                               k=32  e01=-0        (lost: a*a overflowed)
exp(Bivector(NaN e01)).scalar                  1    (the NaN survives in e01)
is_ideal(point with NaN weight)                0
```

**B2.** `g++ -O2`, 2,000,000 chained calls each, run twice:

```text
Cl3 rotate(Vector, Rotor)        :    7.7 ns/call
PGA sandwich(Trivector, Motor)   :  273.4 ns/call
PGA sandwich(Trivector, Rotor)   :  302.0 ns/call
```

**B3.** The same program at `-O0`, at `-O2`, and at
`-O2 -ffp-contract=fast -mfma`, printed with `%a`. `-O0` and `-O2` agree bit
for bit; contraction does not:

```text
-O2          wedge(a, a)  : 0x0p+0 0x0p+0 0x0p+0
-O2 + FMA    wedge(a, a)  : -0x1.e02e5815f3478p-59 0x1.a3b6d3e71d5bp-59 -0x1.52ce06c1f465p-57
-O2          point weight : 0x1p+0
-O2 + FMA    point weight : 0x1.fffffffffffffp-1
```

The FMA build really does contract: `objdump` counts 150 `vfmadd`/`vfmsub`
instructions in the PGA suite's `products.o`. Contraction is also GCC's
default in ISO C++ mode: plain `-std=c++17 -O2 -mfma`, with no contraction
flag, gives the same non-zero `a ^ a`, and adding `-ffp-contract=off` restores
`0x0p+0`.

**C2.** `ga::normalize(r);` and `ga::inverse(r);` as statements compiled with
no diagnostic under `g++ -Wall -Wextra -Wpedantic` and under MSVC `/W4`. They
now fail `tests/compile_fail/`.

**C4.** `constexpr double n = ga::norm(ga::Vector<double>(3, 4, 0)).value;`
is accepted by g++ with `-Wpedantic`. MSVC rejects it:
`error C2131: expression did not evaluate to a constant`, failing at `sqrt` in
`norm.hpp`.

**D3.** The suites were built with
`-O1 -fsanitize=undefined,float-divide-by-zero -fsanitize-undefined-trap-on-error`.
All passed, with no trap. A one-line program calling `inverse(Vector 0)` under
the same flags trapped (`Illegal instruction`), which shows the sanitizer was
live and the original tests simply never divided by zero.
