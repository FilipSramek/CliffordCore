# CliffordCore

CliffordCore is a small, header-only C++17 library for 3D Clifford (geometric)
algebra, Cl(3,0). It provides the full grade ladder — scalars, vectors,
bivectors, trivectors, multivectors and rotors — with the products, involutions
and rotation machinery that go with them.

Everything is `constexpr`, templated over any arithmetic type, and dependency
free.

## Features

**Types**

- `Scalar`, `Vector`, `Bivector`, `Trivector` — the four grades
- `Multivector` — all grades at once
- `Rotor` — a rotation, as a scalar plus a bivector

**Operations**

- Products: dot / left contraction (`|`), wedge (`^`), and the geometric product
  (`*`) — defined for **every** operand pair, backed by the full 8x8
  multiplication table on multivectors
- The rest of the inner-product family: `right_contraction`, `scalar_product`
- Mixed-grade addition and subtraction, plus compound assignment
- `norm`, `squared_norm`, `normalize`
- All three involutions — `reverse`, `involute`, `conjugate` — and `inverse`
- `dual` across all four grades, and `grade0`–`grade3` projection
- Geometry: `reflect`, `project`, `reject`, onto a line or a plane
- Rotations: `exp`, `log`, `sandwich`/`rotate` (for vectors, bivectors,
  trivectors and multivectors), `identity_rotor`, `rotor_from_axis_angle`,
  `rotor_between`, `slerp`
- Exact `==` / `!=`, a separate `approx_equal`, and `operator<<` / `to_string()`
  at round-trip precision

Header-only, `constexpr` throughout except the printing helpers. Take everything
with one umbrella header, or include only the pieces you need — every header
stands alone.

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

Compile from the repository root:

```bash
g++ -std=c++17 -Iinclude your_file.cpp -o your_file
```

## Project structure

- `include/cliffordcore/` — the library, one directory per algebra.
  `cliffordcore/cl3.hpp` pulls in all of Cl(3,0); the six type headers and
  `cl3/operations/` can also be included individually. `cliffordcore/detail/`
  holds the few helpers that are algebra-independent
- `examples/` — five runnable programs, from the basic products through
  rotations to a small geometry toolkit
- `tests/` — the test suite, a single translation unit
- `docs/` — guides, a conventions reference, and a Doxyfile

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

**Version 0.1.0. Pre-1.0: the API is not yet stable.**

CliffordCore is an educational and experimental implementation, suitable for
learning, prototyping, and exploring geometric algebra in C++. It implements
Cl(3,0), 3D Euclidean geometric algebra, and the test suite covers it thoroughly.

**The namespace reorganisation has landed.** Each algebra gets its own
namespace under `CliffordCore`, with the *same* type names in each, so choosing
an algebra is one line and switching is a one-line edit:

| Algebra | Header | Namespace | Status |
| --- | --- | --- | --- |
| Cl(3,0), 3D Euclidean | `<cliffordcore/cl3.hpp>` | `CliffordCore::Cl3` | shipping |
| Cl(3,0,1), 3D projective | `<cliffordcore/pga.hpp>` | `CliffordCore::PGA` | planned |
| Cl(2,0), 2D Euclidean | `<cliffordcore/cl2.hpp>` | `CliffordCore::Cl2` | planned |

That is why the types dropped their dimension suffix — `Vector3` became
`Cl3::Vector`, which said "3" twice, and would have been actively wrong in
`Cl2`. **If you are on an earlier commit, this rename will break your code**;
the fix is mechanical. Further API changes are still possible before 1.0.

**No license yet.** Until a LICENSE file lands the code is all rights reserved,
so this is not yet usable in your own project. That is being resolved.
