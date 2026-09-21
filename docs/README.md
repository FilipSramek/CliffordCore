# CliffordCore documentation

CliffordCore is a header-only C++17 library for Clifford (geometric) algebra.
Each algebra lives in its own namespace under `CliffordCore`, with the **same
type names** in each, so choosing one is a single line:

```cpp
#include <cliffordcore/cl3.hpp>
namespace ga = CliffordCore::Cl3;
```

| Algebra | Header | Namespace | Status |
| --- | --- | --- | --- |
| Cl(3,0), 3D Euclidean | `<cliffordcore/cl3.hpp>` | `CliffordCore::Cl3` | shipping |
| Cl(3,0,1), 3D projective | `<cliffordcore/pga.hpp>` | `CliffordCore::PGA` | shipping |
| Cl(2,0), 2D Euclidean | `<cliffordcore/cl2.hpp>` | `CliffordCore::Cl2` | planned |

**Unless a page says otherwise it documents Cl(3,0).** The type names carry
across to PGA but the meanings do not — a `Vector` is a direction in Cl(3,0) and
a *plane* in Cl(3,0,1) — so PGA has its own conventions page.

## Pages

| Page | What it covers |
| --- | --- |
| [getting-started.md](getting-started.md) | Installing, a first program, choosing `T`, what is `constexpr`, how to print |
| [conventions.md](conventions.md) | **The Cl(3,0) sign rules.** Basis, return types, duality, orientation, known gaps |
| [pga.md](pga.md) | **The Cl(3,0,1) sign rules**, and what a degenerate `e0` changes |
| [rotations.md](rotations.md) | Rotors in practice: building, applying, composing, drift, interpolation |
| [cheatsheet.md](cheatsheet.md) | Return-type matrices and a task-to-call index, on one screen |
| [api-reference.md](api-reference.md) | Every header, every signature |

**If a sign looks wrong, go to [conventions.md](conventions.md) first** — or to
[pga.md](pga.md) if you are in PGA. Almost every surprise in this library traces
back to one of the choices listed there.

## Examples

Six runnable programs in [`../examples/`](../examples/), building in difficulty:

| Example | Shows |
| --- | --- |
| `01_getting_started.cpp` | The three products and what each returns |
| `02_grades_and_multivectors.cpp` | The grade ladder, widening, projection, duality, inverses |
| `03_rotations.cpp` | Rotors, the sign trap, composition, drift |
| `04_interpolation.cpp` | `rotor_between`, `slerp`, recovering axis and angle |
| `05_geometry_toolkit.cpp` | Angles, areas, volumes, projection, reflection |
| `06_pga_rigid_motions.cpp` | Cl(3,0,1): planes, lines and points in space, translators, motors, screws |

The first five are Cl(3,0); the sixth switches algebra and is the one to read
for what PGA buys over it.

Build and run them all:

```bash
./build_examples.sh            # POSIX / Git Bash
.\build_examples.ps1           # PowerShell
```

Or one at a time, with an optional substring filter:

```bash
./build_examples.sh 03         # just the rotations example
./build_examples.sh --no-run   # compile without executing
```

Binaries go to `build/examples/`, which is gitignored. In VS Code, the default
build task compiles whichever example file is open.

## Tests

```bash
./build.sh                     # POSIX / Git Bash
.\build.ps1                    # PowerShell
```

One translation unit per algebra, each its own binary, so the two never meet in
a single TU. A green run prints `911 checks, 0 failed.` for Cl(3,0) and
`1741 checks, 0 failed.` for Cl(3,0,1), then `2 test suites, 0 failed.`
`build.sh` runs `check_headers.sh` first, which compiles all 64 headers alone in
their own translation units.

## Generated API reference

```bash
cd docs && doxygen Doxyfile
```

**The `cd docs` is required** — the Doxyfile's paths are relative to the working
directory, not to the file, so running it from elsewhere silently produces
nothing. Output lands in `docs/api/html/index.html`, and warnings in
`docs/doxygen-warnings.log`. Both are gitignored.

Graphviz is not required; diagram generation is off.
