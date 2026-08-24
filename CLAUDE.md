# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

CliffordCore is a header-only C++17 template library implementing 3D Clifford (geometric) algebra. Everything lives in `namespace CliffordCore`, every type is `template<typename T>` guarded by `static_assert(std::is_arithmetic<T>::value, ...)`, and every operation is `constexpr`. There is no build system, no library artifact, and no third-party dependency — consumers just add `include/` to their include path.

## Build and test

The whole test suite is one translation unit. Any of these compile and run it:

```bash
./build.sh                      # Git Bash / POSIX
.\build.ps1                     # PowerShell (substs the repo root to Z:, builds, runs, unmounts)
g++ -std=c++17 -Iinclude tests/test_core.cpp -o build/test_core && ./build/test_core.exe
```

Toolchain in use is MinGW-W64 g++ 13.2.0 (UCRT) on PATH. `.vscode/tasks.json` hardcodes a *different* compiler (`C:\Strawberry\c\bin\g++.exe`) and builds only the active file — it is not the project build.

Two things to know before running a build:

- **`build/test_core.exe` and `build/test.exe` are tracked in git and `.gitignore` is empty.** Building rewrites `build/test_core.exe`, so `git status` will show it modified afterwards. Restore it with `git checkout -- build/test_core.exe` unless the binary change is intended.
- `build.sh` passes `-o build/test_core`; MinGW appends `.exe`, which is why the tracked filename has the extension.

### Running a single test

There is no test framework. [tests/test_core.cpp](tests/test_core.cpp) is a single `main()` calling local `expect_close` / `expect_true` helpers that print to `stderr` and `std::exit(1)` on the first failure — so a failing run reports one failure, not all of them. To exercise a subset, comment out assertions in `main()` or compile a scratch `.cpp` against `-Iinclude` rather than trying to filter at runtime.

`build/test.cpp` is a leftover scratch file, not part of any target. It includes headers as `"vector3.hpp"` (paths relative to `include/`), so it only compiles with `-Iinclude`.

## Architecture

### The circular-include pattern (most important convention)

`Scalar<T>` multiplies with every other type, and `Vector3<T>` / `Bivector3<T>` multiply back with `Scalar<T>`. That mutual dependency is resolved by a deliberate **two-phase declare-then-define** layout, not by include order luck. Both directions of the cycle follow it:

- [include/scalar.hpp](include/scalar.hpp) — forward-declares `Vector3`, `Bivector3`, `Trivector3`, `Multivector3`, `Rotor3`; declares `Scalar<T>` with cross-type `operator*` as **declarations only**; *then*, after the namespace closes, `#include`s all five headers; then defines those operators out-of-line.
- [include/vector3.hpp](include/vector3.hpp) and [include/bivector3.hpp](include/bivector3.hpp) — mirror image: forward-declare `Scalar<T>`, define same-type ops inline in the class body, declare only the `Scalar`-mixing ops, `#include "scalar.hpp"` after the class, then define them.
- `trivector3.hpp`, `multivector3.hpp`, `rotor3.hpp` don't participate — they include their dependencies normally at the top.

Consequences when editing:

- **Adding an operator that mixes two types requires following this pattern.** Defining it inline in the class body where the other type is only forward-declared produces incomplete-type errors.
- The trailing includes are load-bearing. Do not "tidy" them up to the top of the file.
- Every header is a valid standalone entry point — including just `vector3.hpp` transitively pulls in the scalar operations.
- All headers use `#pragma once`, no include guards.

### Grades and component layout

| Type | Grade | Members |
|---|---|---|
| `Scalar<T>` | 0 | `value` |
| `Vector3<T>` | 1 | `x, y, z` |
| `Bivector3<T>` | 2 | `xy, xz, yz` |
| `Trivector3<T>` | 3 (pseudoscalar) | `e123` |
| `Multivector3<T>` | mixed | `scalar, vector, bivector, trivector` |
| `Rotor3<T>` | 0 + 2 | `scalar, bivector` |

**The bivector basis is `(xy, xz, yz)` — note `xz`, not the more common `zx`.** Every sign convention in the wedge product and dual follows from that choice, so it is the first thing to check when a result's sign looks wrong.

### Operations

Free functions in `include/operations/`, all in `namespace CliffordCore`. Each comes as an **operator plus a named wrapper**, with the operator holding the implementation and the named function delegating to it:

- `a | b` = `dot_product(a, b)` → `Scalar<T>`
- `a ^ b` = `wedge_product(a, b)` → `Bivector3<T>`
- `a * b` (`Vector3` × `Vector3`) = `geometric_product(a, b)` → `Multivector3<T>` (scalar part = dot, bivector part = wedge)
- `rotor_product(a, b)` → `Rotor3<T>`, same components packed as a rotor

Note `operator*` is overloaded across two layers with different meanings: the `Scalar`-mixing forms are **class members** in the type headers, while `Vector3 * Vector3` is a **free function** in [include/operations/geometric_product.hpp](include/operations/geometric_product.hpp). `inverse` is built on `squared_norm`, which is built on the raw components.

### Current coverage

The API is deliberately partial; these are gaps, not bugs to route around:

- `inverse` is defined for `Vector3` only. `norm` / `squared_norm` for `Vector3` and `Bivector3` only. `reverse` and `dual` for `Vector3` and `Bivector3` only.
- The geometric product exists only for `Vector3` × `Vector3`. There is no `Multivector3` × `Multivector3`, no `Rotor3` × `Rotor3`, and no sandwich product / rotor application to a vector.
- `dual()` is a straight component copy (`x→xy, y→xz, z→yz` and back), not multiplication by the pseudoscalar. [tests/test_core.cpp](tests/test_core.cpp) asserts this behavior, so changing it to the algebraic dual means updating those assertions deliberately.
- `examples/` and `docs/` are referenced by the README but are empty.

## Code style

- Doxygen `@brief` / `@param` / `@property` / `@return` on every public type and function. The house style places the doc comment **between** `template<typename T>` and the class or function it documents — match it.
- Constructors and operators are `constexpr`; same-type arithmetic is defined inline in the class body, cross-type arithmetic out-of-line (see the include pattern above).
- Commit messages follow `<Type>, <description>` — e.g. `Feature, add rotor type`, `Bug, removed dot and wedge`, `Docs, added a basic readme file`.
