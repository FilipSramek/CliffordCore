# Coding guidelines

These are the rules for writing code in CliffordCore, and the reason behind each
one. Most of them were learned the hard way: a rule here usually exists because
breaking it once shipped a bug. The library's long-term target is
safety-critical and hard real-time use, so several rules come from
[safety-audit.md](safety-audit.md); the finding each one answers is named in
brackets, like (A1).

**Where a rule and existing code disagree, the code is the thing to fix** —
the audit lists the known cases. A pull request should not add a new one.

## Language and dependencies

- **C++17, header-only, no third-party dependency.** A consumer adds
  `include/` to the include path and nothing else.
- **Standard C++ only.** It must compile cleanly under g++
  `-std=c++17 -Wall -Wextra -Wpedantic` and under MSVC
  `/std:c++17 /W4 /permissive-`. Warnings count as failures even though the
  build does not use `-Werror`.
- **No compiler extensions that change meaning.** One is already in the code:
  g++ will evaluate `std::sqrt` at compile time, MSVC will not, so a
  `constexpr` function that calls `<cmath>` only runs at compile time on GCC
  (C4). Never write a `static_assert` that depends on it.

## The real-time rules

The math path — everything except printing — must contain **no**:

| Not allowed | Why |
| --- | --- |
| heap allocation (`new`, `std::vector`, `std::string`, …) | allocation time is unbounded, and most safety standards forbid it after start-up |
| exceptions (`throw`, `try`) | unwinding time is unbounded, and many targets build with `-fno-exceptions` |
| loops or recursion | each call must run a fixed number of operations, so its worst-case time is obvious from the source |
| virtual functions or RTTI | dispatch costs and `-fno-rtti` builds |
| global or `static` mutable state | it would break reentrancy and thread safety |

Together these make **every call's execution time bounded by construction**,
which is the library's strongest real-time property. One loop or one
allocation breaks that argument for the whole library, so treat these as hard
rules, not preferences.

Printing (`to_string`, `operator<<`, `detail::format_component`) is the only
exception. It allocates, and must stay out of real-time code paths.

## Types

- **Every type is `template<typename T>`, guarded by the same assertion:**

  ```cpp
  static_assert(std::is_floating_point<T>::value,
                "Vector can only be instantiated with floating-point types.");
  ```

  It is `is_floating_point`, not `is_arithmetic`. Integer types used to compile
  and then silently truncate inside `norm`, `inverse` and `exp`. A new type uses
  the same guard and the same message, and gets a check in
  `tests/compile_fail/`.
- **The same type names in every algebra** (`Scalar`, `Vector`, `Bivector`,
  `Trivector`, `Multivector`, `Rotor`), with no dimension suffix, so switching
  algebra is one `namespace ga = ...;` line. An algebra may add types, omit a
  grade it does not have, or add geometric aliases, but never rename.
- **Default constructors zero every member.** Nothing is ever left
  uninitialized.
- **Data members are public**, documented with a trailing `///<` comment.
- **A constructor that loses information is `explicit`**, such as
  `Rotor(const Multivector&)`, which drops grades 1 and 3. While it was
  implicit, `r + m` compiled and silently discarded half the multivector.
  Today several single-argument constructors that do not narrow are still
  implicit (C1), and they cause the `Scalar * T` ambiguity. Do not add another
  implicit one.

## Declaring functions: `[[nodiscard]]`, `constexpr`, `noexcept`

Every function in the math path is declared like this:

```cpp
[[nodiscard]] constexpr Vector operator+(const Vector& other) const noexcept {
    return Vector(x + other.x, y + other.y, z + other.z);
}
```

Each word does a separate job.

### `[[nodiscard]]`

**What it is.** A standard C++17 *attribute*. Attributes are written in double
square brackets before a declaration, like `[[deprecated]]` or
`[[maybe_unused]]`. `[[nodiscard]]` tells the compiler that calling this
function and ignoring its result is almost certainly a mistake, so it should
warn.

**Why it is here.** In this library every operation returns a new value and
leaves its arguments alone. That makes this an easy bug to write:

```cpp
ga::Rotor<double> r = /* ... drifted ... */;
ga::normalize(r);           // looks like it normalizes r; it does nothing
ga::Rotor<double> u = ga::normalize(r);   // what was meant
```

Without the attribute, the second line compiled silently under every warning
flag on both compilers. With it, g++ says *ignoring return value of ...,
declared with attribute 'nodiscard'*, and MSVC reports warning C4834. When
discarding a result really is intended, which it rarely is, say so with
`(void)ga::normalize(r);`.

It changes nothing at run time. It is purely a compile-time check, and maps to
MISRA C++:2023 Rule 0.1.2, "the value returned by a function shall be used".

**Where it goes:** on every function and operator that returns a value. It does
**not** go on:

- constructors, which have no return value;
- functions that return a reference, such as `operator+=` (returns `*this` so
  that `a += b += c` chains) and `operator<<` (returns the stream). Ignoring
  those results is normal.

`tests/compile_fail/*_discarded_*.cpp` turn the warning into an error and
check that it fires, so the attribute cannot quietly disappear.

### `constexpr`

The function can be evaluated at compile time when its arguments are
constants, so `constexpr Vector<double> k(1, 2, 3);` is a compile-time
constant. All arithmetic is `constexpr`. The exceptions are the printing
helpers, `approx_equal`, and the `inline` functions that return a raw `T`
from `std::atan2`. For anything that calls `<cmath>`, compile-time
evaluation works only on GCC (C4).

### `noexcept`

**What it is.** An *exception specification*: a promise that the function never
throws. If something inside it did throw anyway, the exception would not
propagate to the caller; the program would end through `std::terminate()`.

**Why it comes after `const`.** That is the grammar. After the parameter list,
a member function's declarator takes, in this order: `const` (this function
does not modify the object), an optional `&`/`&&`, then `noexcept`. So
`) const noexcept {` is the only legal spelling; `) noexcept const {` is a
syntax error. A free function has no `const` slot, so it reads
`) noexcept {`.

**Why it is here.** Nothing in the math path can throw, since it only does
floating-point arithmetic and calls `<cmath>`. `noexcept` makes that fact part
of the declaration, where a reader, the compiler and a static analyser can all
see it. A safety analysis can then prove there is no exception path at all
(AUTOSAR A15-4-4 asks for exactly this), and a test can check it with
`static_assert(noexcept(ga::normalize(v)))`.

**Where it goes:** on every function and constructor, except the printing ones.
`to_string`, `operator<<` and `format_component` build a `std::string`, which
can throw `std::bad_alloc`; marking them `noexcept` would turn an out-of-memory
condition into an immediate `std::terminate()`.

**One trap:** a member declared in the class and defined outside it must carry
`noexcept` in **both** places, or the build fails. The cross-type operators
follow this pattern, so check both copies.

### The full order

```text
[[nodiscard]] constexpr ReturnType name(parameters) const noexcept
```

- A free function has no `const`.
- A constructor has no return type and no `[[nodiscard]]`.
- A compound assignment returns a reference and has no `[[nodiscard]]`.

## Operators and products

- **`operator*` returns the most general type the product can produce**, and
  there is exactly one `operator*` per pair of operand types. C++ cannot
  overload on return type: a second `operator*(Vector, Vector)` returning
  `Rotor` once made every `v * v` ambiguous.
- **Narrower results are named functions:** `rotor_product`, `rotor_sum`,
  `to_rotor`, `to_multivector`, `motor_product`.
- **Mixed-grade `+` and `-` always widen to `Multivector`**, even when the pair
  happens to be a rotor. Uniform rules keep chained expressions predictable.
- **Every mixed product goes through the general product**,
  `geometric_product(Multivector, Multivector)`, so there is one multiplication
  table per algebra. A faster closed form is allowed, as `Vector * Vector` and
  `Rotor * Rotor` have, but a test must assert that it agrees with the general
  product.
- **Mind the precedence of `^` and `|`.** They bind more loosely than `+`, `-`,
  `*` and `==` (C6):

  | Written | Parsed as |
  | --- | --- |
  | `a ^ b + c` | `a ^ (b + c)` |
  | `a \| b * c` | `a \| (b * c)` |
  | `a ^ b == c` | `a ^ (b == c)` |

  Always parenthesize them, `(a ^ b) + c`. In safety-relevant code, prefer the
  named forms `wedge_product(a, b)` and `left_contraction(a, b)`.

## Header layout

- **`#pragma once`**, then the file's doc block, then the includes.
- **The doc block** carries `@file`, `@brief`, a paragraph of explanation, and
  the house fields. Copy them from a neighbouring header:

  ```cpp
  /**
   * @file cliffordcore/cl3/operations/dual.hpp
   * @brief Duality: multiplication by the pseudoscalar.
   *
   * What the file is for, and the one convention a reader must know.
   *
   * @author Filip Sramek
   * @version 0.2.0
   * @date 2026
   * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
   * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
   */
  ```

- **A header includes everything it uses itself.** In particular, a header
  that calls `std::sqrt` or `std::abs` includes `<cmath>` directly, never by
  relying on another header. That has broken twice. `./check_headers.sh`
  compiles every header alone to catch it.
- **The type headers use a two-phase layout:** declare the class, close the
  namespace, `#include` the headers it mixes with, then define the cross-type
  operators. Those trailing includes are load-bearing, so never move them to
  the top. Adding an operator between two types means following the same
  pattern.
- **Nested namespaces**, `namespace CliffordCore::Cl3 { ... }`.
- **Shared helpers are spelled in full**, `CliffordCore::detail::format_component`.
  Each algebra has its own inner `detail`, which hides the outer one.
- **There is no header that includes every algebra**, and none should be added.
  Mixing two algebras in one expression is a bug, and separate namespaces make
  it a compile error.

## Documentation comments

- **Every public type and function gets `@brief`, `@param`, `@tparam` and
  `@return`**, and `@pre` wherever a precondition exists, such as a unit rotor
  or a non-zero argument.
- **The doc comment sits between `template<typename T>` and the declaration:**

  ```cpp
  template<typename T>
  /**
   * @brief Normalizes a rotor to unit norm.
   * @param r The rotor to normalize.
   * @return The rotor divided by its norm.
   * @pre r is non-zero.
   */
  [[nodiscard]] constexpr Rotor<T> normalize(const Rotor<T>& r) noexcept {
  ```

- **Data members are documented with a trailing `///<`.** Never use
  `@property`: it is not a C++ command, and once silently suppressed every
  class brief in the library. Never use `@license` either; it is not a Doxygen
  command.
- **Doxygen must stay at zero warnings.** Run `cd docs && doxygen Doxyfile`,
  and `docs/doxygen-warnings.log` must be empty. Doxygen also reads the
  Markdown pages in `docs/`. Inside them:
  - do not use in-page links such as `[text](#heading)`, which it cannot
    resolve;
  - write a `|` inside a table cell as `\|`.
- **Say what a degenerate input returns.** If zero, an ideal element or a
  non-unit versor gives a particular result, the `@return` or `@pre` says so.

## Numerics

- **Angles come from `atan2(sin, cos)`, never `acos(cos)`.** `acos` loses
  small angles to rounding; the float property tests caught it.
- **In Cl2, the `+ T(0)` inside `atan2` calls is load-bearing.** It turns `-0.0`
  into `+0.0`, so a half turn reports `+π`, not `-π`.
- **Format numbers through `detail::format_component`**, never
  `std::to_string`, which prints six fixed decimals and shows every 1e-16
  residual in this library as `0.000000`.
- **Compare floating-point values with a tolerance** (`approx_equal`,
  `check_close`). Exact `==` is only for values equal by construction, and
  even then FMA contraction can break it: with `-ffp-contract=fast`, `a ^ a`
  comes out near 1e-17, not 0 (B3). The existing exact `n == T(0)` guards are a
  known weakness (C3); do not add new ones.
- **Respect the input domain:** every component within [2^-30, 2^30] in float
  and [2^-250, 2^250] in double. The bound comes from `a * a` in PGA
  `inverse(Bivector)`, where `a` is already a squared norm. A new expression
  that squares a squared norm narrows the domain, so re-derive it and update
  `domain_max` in the three `special_values.cpp` files and the audit.
- **Do not add a new silent sentinel.** Today a zero axis gives the identity
  rotor and a distance to an ideal point gives `0`. These are known findings
  (A1), to be replaced by one checked policy. Until then, a new function
  documents its degenerate case and pins it in `special_values.cpp`.

## Tests

The full description is in [README.md](README.md) and the test files
themselves. The rules:

- **A test goes in the category file of its kind** (`products.cpp`,
  `rotations.cpp`, …), inside that file's unnamed namespace, and gets one line
  in the file's `run_<category>()`. A test that is never called draws a
  `defined but not used` warning, which is the build's guard against forgetting
  to wire it in.
- **A new public function is called in `instantiation_sweep.cpp`** for float,
  double and long double. Templates only fail when instantiated: `to_string()`
  once shipped broken on five of six types with every test green.
- **A new identity goes in `properties.cpp`**, which runs every type with
  tolerances scaled to its epsilon. Float is where precision bugs show.
- **NaN, zero, ideal and domain-edge behaviour goes in `special_values.cpp`.** A
  new entry point joins the NaN sweep, which asserts that a NaN in gives a
  non-finite result out.
- **Check math against an independent reference** (the general product, the
  blade-product rule, or `std::complex` in Cl2), not against a hand-derived
  number.
- **Never let a test hide a defect.** `project(Line, Plane)` returned the line
  reversed for as long as its test compared `std::abs` of the direction. If a
  result is known to be wrong and not yet fixed, pin the wrong value with a
  comment that names the defect, so the fix has to update it deliberately.

## Before a commit

| Check | Command | When |
| --- | --- | --- |
| Headers standalone, compile-fail checks, all suites | `./build.sh` | every change to code or tests |
| Optimized, release and FMA builds | `./build.sh --variants` | anything numerical |
| Examples | `./build_examples.sh` | API changes |
| Doxygen at zero warnings | `cd docs && doxygen Doxyfile` | header comments or `docs/*.md` |
| MSVC | `cl /std:c++17 /EHsc /W4 /permissive- /bigobj /Iinclude tests\cl3\*.cpp`, per suite | portability-sensitive changes |

Run the MSVC command from a Developer Command Prompt (`vcvars64.bat`). Copy
`include/` and `tests/` to a path without non-ASCII characters first. MSVC
prints its diagnostics in the system language, so search the output for the
code pattern `C####:` rather than for the word "warning".

A green run has 0 failed checks and no warnings, and leaves `git status`
showing only the files you meant to change.

## Commits

- **Small, one concern each.** A mechanical sweep across the library is split
  per algebra.
- **One line, `<Type>, <description>`**, where the type is one of `Feature`,
  `Fix`, `Refactoring`, `Test`, `Docs` or `Build`. For example:
  `Fix, PGA project(Line, Plane) keeps the line's orientation`.
- **If an audit finding is fixed,** update its Status in
  [safety-audit.md](safety-audit.md) in the same commit.
