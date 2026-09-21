# CliffordCore

CliffordCore is a small, header-only C++17 library for 3D Clifford (geometric)
algebra. It ships two algebras, each in its own namespace and using the same
type names, so switching between them is a one-line edit:

- **Cl(3,0)** — 3D Euclidean. Directions and rotations about the origin.
- **Cl(3,0,1)** — 3D projective (PGA). Planes, lines and points *anywhere* in
  space, with translations, rotations about arbitrary axes, and the screw
  motions that combine them.

Everything is `constexpr`, templated over any arithmetic type, and dependency
free.

## Features

**Cl(3,0) types** — `namespace CliffordCore::Cl3`

- `Scalar`, `Vector`, `Bivector`, `Trivector` — the four grades
- `Multivector` — all grades at once
- `Rotor` — a rotation, as a scalar plus a bivector

**Cl(3,0,1) types** — `namespace CliffordCore::PGA`

- `Scalar`, `Vector`, `Bivector`, `Trivector`, `Quadvector` — the five grades,
  which here mean a number, a **plane**, a **line**, a **point** and the
  pseudoscalar. `Plane`, `Line` and `Point` are aliases that say so
- `Multivector` — all sixteen components
- `Rotor` — a rotation about an axis through the origin
- `Translator` — a translation
- `Motor` — any rigid motion: a screw. `Twist` aliases `Bivector`, which is
  what `exp` takes and `log` returns

**Operations, in both algebras**

- Products: dot / left contraction (`|`), wedge (`^`), and the geometric product
  (`*`) — defined for **every** operand pair, backed by one full multiplication
  table on multivectors (8x8 in Cl3, 16x16 in PGA)
- The rest of the inner-product family: `right_contraction`, `scalar_product`
- Mixed-grade addition and subtraction, plus compound assignment
- `norm`, `squared_norm`, `normalize`
- All three involutions — `reverse`, `involute`, `conjugate` — and `inverse`
- `dual`, and grade projection across every grade
- `exp`, `log`, `sandwich`, `slerp`
- Exact `==` / `!=`, a separate `approx_equal`, and `operator<<` / `to_string()`
  at round-trip precision

**Cl(3,0) adds** `reflect`, `project`, `reject` onto a line or a plane;
`identity_rotor`, `rotor_from_axis_angle`, `rotor_between`.

**Cl(3,0,1) adds** the `meet` and `join` of planes, lines and points;
factories — `plane`, `point`, `line_through_points`, `line_from_planes` — and
the accessors that read them back; `translator`, `motor_from_line_angle`,
`screw`, `motor_between`; `reflect` and `project` for objects anywhere in space;
and `distance` and `angle` for every pair of plane, line and point.

Header-only, `constexpr` throughout except the printing helpers. Take one
algebra with one umbrella header, or include only the pieces you need — every
header stands alone. There is deliberately no header that pulls in both.

## Usage

```cpp
#include <iostream>
#include <cliffordcore/cl3.hpp>

namespace ga = CliffordCore::Cl3;

int main()
{
    ga::Vector<double> a(1, 2, 3);
    ga::Vector<double> b(4, 5, 6);

    // The geometric product keeps both the dot and the wedge.
    auto product = a * b;
    std::cout << "dot:   " << product.scalar.value << "\n";
    std::cout << "wedge: " << product.bivector.xy << ", "
              << product.bivector.xz << ", " << product.bivector.yz << "\n";

    // A quarter turn about +z, applied to the x axis. Rotors compose by
    // multiplication and interpolate with slerp.
    auto turn = ga::rotor_from_axis_angle(
        ga::Vector<double>(0, 0, 1), 1.5707963);
    auto spun = ga::rotate(ga::Vector<double>(1, 0, 0), turn);
    std::cout << "e1 turned: (" << spun.x << ", " << spun.y << ", " << spun.z << ")\n";
}
```

The same program in PGA is one changed include and one changed alias, but the
objects mean different things — and can sit anywhere in space:

```cpp
#include <cliffordcore/pga.hpp>

namespace ga = CliffordCore::PGA;

int main()
{
    // A rotation about a line that does NOT pass through the origin --
    // something a Cl(3,0) rotor cannot express at all.
    auto axis = ga::line_through_points(ga::point(1.0, 0.0, 0.0),
                                        ga::point(1.0, 0.0, 1.0));
    auto turn = ga::motor_from_line_angle(axis, 3.1415927);
    auto moved = ga::sandwich(ga::point(2.0, 0.0, 0.0), turn);   // the origin

    // One sandwich moves planes, lines and points alike.
    auto lifted = ga::sandwich(ga::plane(0.0, 0.0, 1.0, 0.0),
                               ga::translator(0.0, 0.0, 5.0));   // z = 5

    auto p = ga::position(moved);
    std::cout << "(2,0,0) swung round to (" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";
    std::cout << "ground lifted to z = " << -lifted.e0 << "\n";
}
```

Compile from the repository root:

```bash
g++ -std=c++17 -Iinclude your_file.cpp -o your_file
```

## Project structure

- `include/cliffordcore/` — the library, one directory per algebra.
  `cliffordcore/cl3.hpp` pulls in all of Cl(3,0) and `cliffordcore/pga.hpp` all
  of Cl(3,0,1); the type headers and each `operations/` directory can also be
  included individually. `cliffordcore/detail/` holds the few helpers that are
  algebra-independent
- `examples/` — six runnable programs, from the basic products through
  rotations to a small geometry toolkit, and one on PGA rigid motions
- `tests/` — one translation unit per algebra; `./build.sh` builds and runs
  every `tests/*.cpp`
- `docs/` — guides, a conventions reference per algebra, and a Doxyfile

## Building

```bash
./build.sh                      # run the tests
./build_examples.sh             # build and run every example
cd docs && doxygen Doxyfile     # generate the API reference
```

PowerShell equivalents: `.\build.ps1` and `.\build_examples.ps1`.

`./build_examples.sh 03` builds just the matching example; `--no-run` compiles
without executing.

## Documentation

Start with [docs/README.md](docs/README.md). The page worth reading before
anything else is [docs/conventions.md](docs/conventions.md) — the bivector basis
is `(xy, xz, yz)`, the dual negates when applied twice, and `exp(theta*B)`
rotates by `2*theta`. Nearly every surprise traces back to one of those.

## Status

**Version 0.2.0. Pre-1.0: the API is not yet stable.**

CliffordCore is an educational and experimental implementation, suitable for
learning, prototyping, and exploring geometric algebra in C++. Each algebra gets
its own namespace under `CliffordCore`, with the *same* type names in each, so
choosing an algebra is one line and switching is a one-line edit:

| Algebra | Header | Namespace | Status |
| --- | --- | --- | --- |
| Cl(3,0), 3D Euclidean | `<cliffordcore/cl3.hpp>` | `CliffordCore::Cl3` | shipping |
| Cl(3,0,1), 3D projective | `<cliffordcore/pga.hpp>` | `CliffordCore::PGA` | shipping |
| Cl(2,0), 2D Euclidean | `<cliffordcore/cl2.hpp>` | `CliffordCore::Cl2` | planned |

That is why the types dropped their dimension suffix — `Vector3` became
`Cl3::Vector`, which said "3" twice, and would have been actively wrong in
`Cl2`. Both algebras are covered thoroughly by their own test suite, which
checks the multiplication table against an independent implementation of the
basis-blade rule rather than against memorised numbers.

The names carry across but the *meanings* do not: a `Vector` is a direction in
Cl(3,0) and a plane in PGA, and none of Cl(3,0)'s sign conventions survive the
degenerate metric. [docs/conventions.md](docs/conventions.md) and
[docs/pga.md](docs/pga.md) are separate pages for exactly that reason. Further
API changes are still possible before 1.0.

## License

Apache License 2.0 -- see [LICENSE](LICENSE) for the full text. You may use,
modify and redistribute the library, including commercially, provided you keep
the notice and state your changes; the patent grant and the "AS IS" disclaimer
come with it. Every header carries `SPDX-License-Identifier: Apache-2.0`.
