# Aerospace certification audit: DO-178C and MISRA C++:2023

**Date:** 2026-10-21 · **Audited:** commit `fb13578` · **Companion to:**
[safety-audit.md](safety-audit.md), whose finding IDs (A1, C4, …) are reused
here where the same problem appears.

This audit asks one question: what would stop CliffordCore from being used in
airborne software certified under **DO-178C**, with **MISRA C++:2023** as the
coding standard? It covers the code (what a MISRA checker reports) and the
evidence a certification would require (what DO-178C and its supplements
expect to exist).

Unlike the first audit, the MISRA results here come from a **real rule
checker**, not from reading the rules:

| Tool | Used for |
| --- | --- |
| Polyspace Bug Finder R2024b, `-misra-cpp-2023 all` plus default defect checkers | every MISRA C++:2023 rule Polyspace implements, on all three algebras |
| Polyspace Code Prover R2024b | proving, or failing to prove, the absence of run-time errors for full-range inputs |
| gcov (GCC 13.2) | statement and branch coverage of the headers by the existing test suites |
| g++ with `-Wconversion -Wfloat-equal -Wold-style-cast -Wshadow -Wuseless-cast -Warith-conversion` and more | a second opinion on conversions and comparisons |

Polyspace analysed each algebra's `tests/<algebra>/instantiation_sweep.cpp`,
which calls every public function for float, double and long double, so every
template body was instantiated and checked. The commands are under
*Reproducing the results* at the end.

## Verdict

**CliffordCore is not certifiable today, but nothing in it is unfixable.**

- **The code side is moderate work.** Polyspace reports 1,351 MISRA C++:2023
  findings in the library, but they come from about a dozen causes. Most are
  mechanical: replace `T(0)` with `static_cast<T>(0)`, replace `#pragma once`
  with include guards, ref-qualify compound assignments. Three need design
  decisions: the asserts, the degenerate-input policy (safety-audit A1), and
  the implicit constructors (C1). Polyspace found **no** dynamic memory, no
  unused or uninitialized variables, no shadowing, no unbraced `if`, and no
  integer-type misuse in the library.
- **The evidence side is the large part.** DO-178C needs plans, standards,
  requirements, traceability, requirements-based tests, structural coverage
  (up to MC/DC at DAL A), DO-332 analyses for templates and overloading, tool
  qualification and configuration management. Almost none of it exists yet.
  That is normal for a library that was not developed under a certification
  process, and DO-178C has a route for it: producing the lifecycle data
  afterwards for "previously developed software".

## What already helps

These properties make the remaining work smaller than usual, and are worth
protecting:

- **No dynamic memory, exceptions, inheritance, virtual functions or global
  state** in the math. DO-332 (the object-oriented supplement) then needs no
  analysis for inheritance, dynamic dispatch, exception management or memory
  management, which are its most expensive areas.
- **Data and control coupling are trivial.** Every function communicates only
  through parameters and return values, so the coupling analysis DO-178C asks
  for (6.4.4.2.c) is short.
- **Statement coverage is already 99%.** Most missing coverage is in one place
  (the asserts) and has one cause.
- **Header-only templates leave no dead code in the binary.** A function the
  airborne program never calls is never instantiated, so it generates no
  object code and needs no coverage there.
- **Independent oracles in the tests:** the Cayley tables are checked against
  the blade-product rule, and Cl2 against `std::complex`. Those are good
  verification techniques to cite.

## Part 1: MISRA C++:2023

### Results

Library findings only, deduplicated across the three runs. The test code has
228 more, listed separately below. *Required* rules must be met, or a formal
deviation raised for each violation. *Advisory* rules can be deviated from,
or disapplied in a recorded re-categorization plan.

| # | Rule | Category | Library findings | Cause | Fix |
| --- | --- | --- | --- | --- | --- |
| M1 | 8.2.2 C-style and functional-notation casts | Required | 173 | `T(0)`, `T(2)`, `T(0.9)` everywhere a literal meets `T` | `static_cast<T>(0)`, or a small named constant: `zero<T>()` / `one<T>()` |
| M2 | 19.2.1 prevent multiple inclusion | Required | 87 | `#pragma once` is not standard C++, so it does not count | include guards, `#ifndef CLIFFORDCORE_CL3_VECTOR_HPP` |
| M3 | 7.0.2 no conversion to bool | Required | 19 | the assert idiom `assert(n != T(0) && "message")`: the string literal converts to `bool` | goes away with the assert (M10) |
| M4 | 8.2.3 cast removing const | Required | 5 | `(void)q;` used to silence an unused parameter in PGA `norm.hpp` and `quadvector.hpp` | leave the parameter unnamed, or `[[maybe_unused]]` |
| M5 | 15.1.3 single-argument constructors explicit | Required | 3 (13 by inspection) | `Scalar(T)` is flagged. By the rule's text the other ten are too: `Trivector(const T&)`, `Quadvector(const T&)`, Cl2 `Bivector(const T&)`, and the widening `Multivector(Rotor)`, `Multivector(Translator)`, `Multivector(Motor)`, `Motor(Rotor)`, `Motor(Translator)` | safety-audit C1; `explicit`, plus named widening functions where implicit widening is wanted |
| M6 | 0.1.2 return value shall be used | Required | 1 | `stream.precision(...)` in `detail/format.hpp` discards the previous precision | `(void)` is not allowed for this; restructure, or deviate (printing is not flight code) |
| M7 | 7.0.6 numeric assignment | Required | 1 | `max_digits10` (`int`) to `std::streamsize` in `format.hpp` | `static_cast<std::streamsize>(...)` |
| M8 | 8.0.1 parentheses | Advisory | 666 | `a * b + c * d` and similar, mostly the Cayley tables (342 in `geometric_product.hpp`) and the sandwich closed forms | parenthesize, or record a re-categorization: the Cayley tables are generated and checked against the blade rule |
| M9 | 16.6.1 symmetric operators as non-members | Advisory | 108 | `operator+`, `operator-`, `operator*`, `operator/` written as members inside the classes | move to free functions next to the cross-type ones |
| M10 | 19.6.1 no `#pragma` | Advisory | 87 | `#pragma once` | include guards (same change as M2) |
| M11 | 6.8.4 ref-qualify members returning `*this` | Advisory | 60 | `operator+=`, `-=`, `*=`, `/=` | `Vector& operator+=(const Vector&) & noexcept` |
| M12 | 8.14.1 no side effects on the right of `&&`/`\|\|` | Advisory | 40 | the right operand calls a function: `approx_equal(a.vector, …)`, `norm(b).value == T(0)` | evaluate into named `const bool`s first, then combine |
| M13 | 19.0.3 `#include` only after directives or comments | Advisory | 26 | the two-phase declare-then-define layout's trailing includes | documented deviation: the layout is load-bearing (see `CLAUDE.md`) |
| M14 | 10.0.1 one variable per declaration | Advisory | 24 | `const T a0 = …, a1 = …;` in the three `geometric_product.hpp` | one per line |
| M15 | 8.2.7 cast pointer to integral | Advisory | 19 | inside MinGW's `assert` macro expansion | goes away with the assert |
| M16 | Dir 0.3.1 floating-point used appropriately | Advisory | 16 | Polyspace's *absorption* defect in PGA `Motor * Motor`: for the sweep's inputs, some float terms are too small to change the sum | not a bug; requires stated accuracy requirements (D3) |

**Not reported by Polyspace for the library:** 6.9.2 (raw integer types), 21.6.1
(dynamic memory), 9.3.1 (unbraced bodies), 11.3.1 (C arrays), 11.6.1
(uninitialized variables), 0.1.1 and 0.2.x (unused values and parameters),
6.4.1 (shadowing), 9.5.x (loops). The test code violates several of these;
the library violates none.

### Rules Polyspace cannot decide, found by review

A checker reports what it can prove from the source. These need a person.

| # | Rule | Finding |
| --- | --- | --- |
| M17 | 4.1.3 no undefined behaviour (Required) | **ISO C++ makes floating-point division by zero, and any overflow, undefined behaviour** ([expr.mul]/4, [expr.pre]/4). The library relies on IEEE 754 instead: Cl3 and Cl2 `inverse(0)` divides by zero on purpose (safety-audit A3), and inputs outside the domain overflow (A4). An implementation that follows IEEE 754 defines the result, but nothing in the library asserts that. Add `static_assert(std::numeric_limits<T>::is_iec559)` next to the floating-point guard, record the IEEE 754 dependency as a deviation, and remove the deliberate divisions by zero. |
| M18 | 4.1.1 conform to ISO C++17 (Required) | **At least 36 `constexpr` functions call `<cmath>` unconditionally** (`norm`, `normalize`, `exp`, …). Under C++17 a `constexpr` function that can never be a constant expression is ill-formed, no diagnostic required. GCC accepts it; MSVC rejects compile-time use (safety-audit C4). Drop `constexpr` from those functions. |
| M19 | Dir 0.3.2 respect function preconditions (Required) | `std::atan2(0, 0)` may raise a domain error by the C standard; Cl2 `angle(zero)`, `angle_between(zero, v)` and `rotor_angle(zero)` call it. `std::sin`/`std::cos` of ±∞ is a domain error, reached by `exp(infinite bivector)`. Each call needs a guard, or a documented precondition the caller satisfies. |
| M20 | 18.5.2 no program-terminating functions (Advisory) | `assert` calls `std::abort()` in 19 places (every `normalize`). Terminating is rarely an acceptable airborne failure mode, and it creates a debug/release divergence (safety-audit A2). Replace it with the checked-result policy from safety-audit Phase 2. |
| M21 | Dir 0.3.1 floating-point used appropriately (Advisory) | **106 source lines compare floating-point values with `==`** (g++ `-Wfloat-equal`), in guards such as `n == T(0)` and in the public `operator==`. Each needs a justification or a tolerance; safety-audit C3 and A8 show what the exact guards miss. |

### The MISRA compliance process

Conformance to MISRA C++:2023 is claimed through the **MISRA Compliance:2020**
framework, and none of its four artefacts exists yet:

1. **A guideline enforcement plan:** which tool checks which rule, and how
   the rules no tool can check are reviewed. Polyspace covers most of them;
   M17–M21 are manual.
2. **A re-categorization plan:** for example, disapplying 8.0.1 for the
   generated Cayley tables, or raising an advisory rule to required.
3. **Deviation records:** one per accepted violation of a required rule, with
   its justification. Likely candidates are M13 (two-phase includes), M17
   (IEEE 754) and M6 (printing).
4. **A compliance summary** for the delivered baseline.

**Test code:** DO-178C does not require test code to follow the coding
standard, but test procedures must be reviewed, so a clean harness helps. Its
228 findings are mainly casts (75), raw `int` (36) and `std::string` (32).

## Part 2: DO-178C

### Structural coverage, measured

Measured with gcov over the existing suites (default build: `-O0`, asserts
on). A line or branch inside a template counts as covered if any
instantiation (float, double or long double) covered it. The *per
instantiation* column counts each instantiation separately, which is the
stricter reading.

| Algebra | Statements | Branch outcomes | Branch outcomes, per instantiation |
| --- | --- | --- | --- |
| Cl2 | 608/613 (99.2%) | 80/99 (80.8%) | 216/297 (72.7%) |
| Cl3 | 808/816 (99.0%) | 113/147 (76.9%) | 321/441 (72.8%) |
| PGA | 1574/1587 (99.2%) | 271/352 (77.0%) | 735/1056 (69.6%) |

What DO-178C requires (Table A-7):

| DAL | Failure condition | Coverage required |
| --- | --- | --- |
| A | catastrophic | statement, decision **and MC/DC**, plus source-to-object traceability |
| B | hazardous | statement and decision |
| C | major | statement |
| D, E | minor, none | none |

Statement coverage nearly meets DAL C today. Decision coverage, needed for B,
falls short for four distinct reasons:

1. **The asserts.** Each `assert` is a decision whose "fail" outcome aborts the
   program, so no in-process test can cover it. Every `normalize` in every
   algebra is affected. Coverage would need either death tests (one process
   per case) or removing the asserts. Removing them is better (M20).
2. **The release-only branch.** `if (n == T(0)) return v;` in `normalize`
   runs only under `NDEBUG`, so the debug build leaves it uncovered. The flight
   configuration is the one that must be measured, so coverage has to be taken
   from a `-DNDEBUG` build, which `./build.sh --variants` already builds.
3. **Short-circuit outcomes in comparisons.** `a.x == b.x && a.y == b.y && …`
   is never tested with *only* `y` (or only `z`) differing. All of
   `comparison.hpp` and `is_ideal` show this. It is also exactly the gap MC/DC
   exposes.
4. **Behaviour that has no test at all:**
   - Cl3 `slerp` between rotors pointing into opposite hemispheres. The
     "negate and take the short way round" branch never runs.
   - Cl3 and Cl2 `rotor_between` of exactly opposite vectors whose `x`
     component is ≥ 0.9. The second fallback axis is never chosen.
   - PGA `project(point, ideal line)`, `inverse` of a zero `Motor`,
     `translator_between` with an ideal point, `motor_from_line_angle` with a
     zero axis, and the zero-norm branches of `angle(plane, plane)` and
     `angle(line, line)`.

**MC/DC.** 99 lines in the headers hold decisions with two or more
conditions (`&&` or `||`), each needing MC/DC test pairs at DAL A. GCC 13
cannot measure MC/DC. GCC 14 adds `-fcondition-coverage`. Certification
projects normally use a qualified tool such as VectorCAST, LDRA or RapiCover.

**Templates and coverage.** Coverage has to be measured on the instantiations
the airborne program actually uses, in the target build. A float flight
program needs the float instantiations covered; the long double tests earn it
nothing. At DAL A, inlining and template instantiation also make
source-to-object traceability (DO-178C 6.4.4.2.b) a real analysis, because
the compiler generates code that has no one-to-one source line.

### Findings

Severity: **B** = blocks certification · **M** = major work · **m** = minor.

| # | Sev | DO-178C / DO-332 / DO-330 | Finding |
| --- | --- | --- | --- |
| D1 | B | Planning (§4, Table A-1) | **No plans or standards.** There is no PSAC, SDP, SVP, SCMP or SQAP, and no requirements, design or code standards. `docs/coding-guidelines.md` is a good start on a Software Code Standard, and should name MISRA C++:2023 and its deviations. |
| D2 | B | Requirements (§5.1, §5.2) | **No high-level or low-level requirements.** The docs describe conventions and behaviour, but nothing is written as verifiable, uniquely identified requirements. The degenerate-input sentinels (safety-audit A1) would become *derived requirements*, which go to the system safety assessment, and an identity rotor for a zero axis is unlikely to pass it. |
| D3 | B | Accuracy requirements | **No stated accuracy.** The tests use tolerances like 64 ε chosen per test. A requirement must state the bound, for example "`rotate` returns each component within N ulp of the exact result for inputs in the domain", and tests must verify that bound, including at the domain edges. Polyspace's absorption finding (M16) is the kind of effect such a bound must account for. |
| D4 | B | Verification (§6.4, Table A-7) | **Tests are not requirements-based or traced.** Structural coverage is only accepted as a *measure* of requirements-based testing, so each test must trace to a requirement. The suites are well structured for this. They need IDs, and a trace matrix. |
| D5 | B | Structural coverage (§6.4.4.2) | Decision coverage 77–81% (B fails), MC/DC unmeasured (A fails). See above. |
| D6 | M | Robustness testing (§6.4.2.2) | `special_values.cpp` covers NaN, zero and domain edges, but robustness cases must be derived from the requirements: out-of-range inputs, invalid inputs, and the boundaries of every equivalence class. Today's tests pin current behaviour rather than required behaviour. |
| D7 | M | DO-332 type conversion | **13 implicit single-argument constructors** (M5). DO-332 asks for an analysis showing that implicit conversions cannot lose information or select an unintended overload. Making them `explicit` removes most of that analysis. |
| D8 | M | DO-332 overloading | **Heavy overloading.** There are 200 `operator*` overloads and up to 20 overloads per function name (`norm`, `reverse`, `normalize`, …). DO-332 requires showing that overload resolution always selects the intended function. The existing `static_assert`s on the return types of all 25, 64 and 81 operand pairs are good evidence to build that argument on. |
| D9 | M | DO-332 parametric polymorphism | **Templates.** Each instantiation used in flight must be verified, not just the template. Restrict the certified configuration to the one type the target uses (normally `float` or `double`) and verify that one. |
| D10 | M | Deactivated code (§6.4.4.3) | `to_string`, `operator<<`, `format_component` and `approx_equal`, plus the debug-only assert branch. In a template library, uninstantiated code is absent from the binary, but the flight build must be *shown* not to instantiate the printing functions. Moving printing out of the type headers (safety-audit B1) turns that into a simple include check. |
| D11 | M | Standard library (§11, §12) | **The C++ standard library is part of the airborne software.** `<cmath>` (`sqrt`, `sin`, `cos`, `atan2`, `abs`) runs on the target and needs a certified math library, from the RTOS or compiler vendor, with lifecycle data. `<string>`, `<sstream>` and `<ostream>` must be absent from the flight build. `<type_traits>`, `<limits>` and `<array>` are compile-time or trivial, but still belong to the vendor library being certified. |
| D12 | M | Floating point on the target | **Floating-point exceptions may trap.** Many avionics processors (PowerPC with FP exceptions enabled, some Cortex-R configurations) trap on division by zero or invalid operations instead of producing Inf or NaN. On such a target, Cl3 `inverse(0)` stops the program. The library must either never raise an FP exception for any input in its domain, or document the required FPU configuration. FMA contraction (safety-audit B3), `FLT_EVAL_METHOD` and double-precision emulation on single-precision FPUs belong in the same configuration statement. |
| D13 | M | Tool qualification (DO-330) | **Every tool whose output is trusted without review needs qualification.** The coverage tool and the static analyser are TQL-5; Polyspace offers a DO-178C qualification kit. `build.sh` and the test harness count as verification tools if their pass/fail verdict is accepted without review. The compiler is not qualified, but its output must be verified, and the exact compiler version and options become configuration items. |
| D14 | M | Target compiler | **Development uses one host compiler.** Certification credit is for the target compiler and options. C++17 is required (97 nested-namespace definitions, 962 `[[nodiscard]]`), and qualified avionics toolchains vary in C++17 support. Check the chosen one early; if it only has C++14, those two features are mechanical to remove. |
| D15 | M | Configuration management (§7, Table A-8) | **No baselines.** There are no tags, no problem-reporting process, and the version disagrees: 27 Cl2 headers say 0.3.0, 64 headers and the version macros say 0.2.0 (safety-audit D5). |
| D16 | M | Independence (Table A-6, A-7) | **One developer.** At DAL A and B, reviews and verification of several objectives must be done by someone other than the author. |
| D17 | m | Previously developed software (§12.1) | The library was developed outside a DO-178C process. That is allowed: §12.1 covers producing the missing lifecycle data afterwards. Every finding here is part of that data. The Apache-2.0 licence places no restriction on certified use. |

### Run-time error proof (Polyspace Code Prover)

Code Prover tries to *prove*, for every possible input, that no division by
zero, overflow or other run-time error can happen. Each check comes out
green (proven safe), red (proven to fail), orange (not proven either way) or
gray (unreachable). It ran on a driver that feeds unknown double values into
`inverse`, `normalize`, `rotor_from_axis_angle`, `rotate`, `rotor_between`,
PGA `point`, `plane`, `distance` and `inverse(Bivector)`:

| Inputs | Red | Orange | Gray | Green |
| --- | --- | --- | --- | --- |
| any double, including Inf and NaN | 0 | 152 (147 overflow, 2 division by zero, 3 assert) | 3 | 1,671 |
| limited to \|x\| ≤ 2^250 | 0 | 95 (89 overflow, 2 division by zero, 3 assert) | 3 | 1,742 |

**No red: no input was proven to fail.** But the oranges do not go away, and
the reason is the main lesson for DO-178C robustness:

- **Code Prover reasons with intervals.** It can express |x| ≤ 2^250, but not
  "zero or at least 2^-250", a range with a hole in it, and not "this rotor has
  unit length". The remaining overflow oranges are divisions by norms that
  could be arbitrarily small: `normalize`, `inverse`, and the rotor that
  `rotor_from_axis_angle` builds and `rotate` then applies. A proof needs the
  library to **check its preconditions in code**, so the analyser can follow
  the branch. That is the checked API safety-audit Phase 2 recommends; without
  it, each orange must be justified by hand.
- **The two division-by-zero oranges are real.**
  - Cl3 `inverse(Vector)` divides by `squared_norm(v)` with no guard
    (M17, safety-audit A3).
  - PGA `inverse(Bivector)` guards `a == 0`, then divides by `a * a`. For a
    small non-zero line, `a * a` underflows to exactly zero, so the guard does
    not protect the division. Measured: a line of direction 1e-100 gives an
    inverse with `e01 = -inf`, 1e-160 gives NaN, and UBSan's
    `float-divide-by-zero` traps. This is just outside the documented domain
    (2^-250 ≈ 5.5e-76), so it is consistent with safety-audit A4. But the guard
    *looks* protective and is not: on a target with FP traps enabled (D12),
    this input stops the program.
- **The 3 asserts are orange,** because Code Prover cannot show `n != 0` from
  the inputs.
- **The 3 gray blocks** are `normalize`'s release-only branch, which is
  unreachable while asserts are enabled. In DO-178C terms, that is dead code in
  the debug configuration and live code in the release one (D10, M20).

Code Prover is the usual tool for this kind of evidence in DO-178C projects,
and MathWorks offers a qualification kit for it. Once the checked API exists,
re-running this driver is how its run-time safety gets demonstrated.

## Recommended order of work

Each step is independent of the target standard's final choice of DAL, and
each one shrinks the next.

1. **Decide the error-handling policy** (safety-audit Phase 2, item 1). It
   removes the asserts (M3, M15, M20 and the assert coverage gap), the
   deliberate divisions by zero (M17, D12) and the sentinels that would become
   derived requirements (D2). This is the one design decision everything else
   waits on.
2. **Mechanical MISRA fixes**, one small commit per rule per algebra: casts
   (M1), include guards (M2, M10), unused-parameter casts (M4), ref-qualifiers
   (M11), one declaration per line (M14), named conditions (M12), `explicit`
   (M5), and `constexpr` on `<cmath>` callers (M18). Add
   `static_assert(is_iec559)` (M17), and make PGA `inverse(Bivector)`'s guard
   test the value it actually divides by, `a * a`.
3. **Close the test gaps** found by coverage: `slerp` across hemispheres,
   `rotor_between`'s second fallback axis, the untested PGA degenerate cases,
   and each condition of the comparison operators.
4. **Write the MISRA compliance artefacts** (enforcement plan,
   re-categorization plan, deviation records), and add the Polyspace run to
   the build so new violations are caught.
5. **Write requirements** for one algebra first, probably Cl3 or PGA, with
   stated accuracy bounds, and trace the existing tests to them.
6. **The planning documents and tool qualification**, once a program and a
   DAL are known. These belong to the certification project, not the library.

## Reproducing the results

The source tree was copied to an ASCII-only path first, because the repository
path contains `Š`.

```bash
PS="/c/Program Files/Polyspace/R2024b/polyspace/bin"
# MISRA C++:2023, findings in the library headers included:
"$PS/polyspace-bug-finder.exe" -sources tests/cl3/instantiation_sweep.cpp \
    -I include -I tests/cl3 -lang CPP -cpp-version cpp17 -compiler gnu12.x \
    -target x86_64 -misra-cpp-2023 all -checkers default \
    -generate-results-for custom=<absolute path to include> \
    -results-dir res_cl3
"$PS/polyspace-results-export.exe" -format csv -results-dir res_cl3 -output-name res_cl3.csv
```

Without `-generate-results-for`, Polyspace silently reports nothing for headers
found through `-I`. The first run of this audit showed 0 library findings for
that reason.

```bash
# Code Prover, on a driver whose main() reads unknown values from a
# volatile array and calls the functions under test:
"$PS/polyspace-code-prover.exe" -sources cp_driver.cpp -I include -lang CPP \
    -cpp-version cpp17 -compiler gnu12.x -target x86_64 \
    -generate-results-for all-headers -results-dir res_cp
```

```bash
# Coverage: build a suite with --coverage, run it, then read the JSON per object.
for s in tests/cl3/*.cpp; do
    g++ -std=c++17 -O0 --coverage -fprofile-abs-path -Iinclude -c "$s" -o obj/$(basename "$s" .cpp).o
done
g++ --coverage obj/*.o -o test_cl3 && ./test_cl3
(cd obj && for o in *.o; do gcov -t -j -b "$o"; done)   # then merge the JSON per header line
```

```bash
# Floating-point equality and conversions:
g++ -std=c++17 -fsyntax-only -Iinclude -Wall -Wextra -Wpedantic -Wconversion \
    -Wsign-conversion -Wdouble-promotion -Wfloat-equal -Wold-style-cast -Wshadow \
    -Wuseless-cast -Warith-conversion tests/cl3/instantiation_sweep.cpp
```

The g++ run reports only `-Wfloat-equal` in the library: no conversion,
old-style-cast or shadowing warnings.
