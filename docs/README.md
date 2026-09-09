# CliffordCore documentation

CliffordCore is a header-only C++17 library for 3D Clifford (geometric) algebra,
Cl(3,0).

## Pages

| Page | What it covers |
| --- | --- |
| [getting-started.md](getting-started.md) | Installing, a first program, choosing `T`, what is `constexpr`, how to print |
| [conventions.md](conventions.md) | **The sign rules.** Basis, return types, duality, orientation, known gaps |
| [rotations.md](rotations.md) | Rotors in practice: building, applying, composing, drift, interpolation |
| [cheatsheet.md](cheatsheet.md) | Return-type matrices and a task-to-call index, on one screen |
| [api-reference.md](api-reference.md) | Every header, every signature |

**If a sign looks wrong, go to [conventions.md](conventions.md) first.** Almost
every surprise in this library traces back to one of the choices listed there.

## Examples

Five runnable programs in [`../examples/`](../examples/), building in difficulty:

| Example | Shows |
| --- | --- |
| `01_getting_started.cpp` | The three products and what each returns |
| `02_grades_and_multivectors.cpp` | The grade ladder, widening, projection, duality, inverses |
| `03_rotations.cpp` | Rotors, the sign trap, composition, drift |
| `04_interpolation.cpp` | `rotor_between`, `slerp`, recovering axis and angle |
| `05_geometry_toolkit.cpp` | Angles, areas, volumes, projection, reflection |

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

A green run prints `779 checks, 0 failed.`

## Generated API reference

```bash
cd docs && doxygen Doxyfile
```

**The `cd docs` is required** — the Doxyfile's paths are relative to the working
directory, not to the file, so running it from elsewhere silently produces
nothing. Output lands in `docs/api/html/index.html`, and warnings in
`docs/doxygen-warnings.log`. Both are gitignored.

Graphviz is not required; diagram generation is off.
