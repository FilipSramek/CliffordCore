# Getting started

## Requirements

Any C++17 compiler. There is no build system, no library to link, and no
third-party dependency. Development happens with MinGW-W64 g++ 13.2.0; the test
suite is clean under `-Wall -Wextra -Wpedantic` on both C++17 and C++20.

## Integration

Add `include/` to your include path. Then either take everything:

```cpp
#include <cliffordcore.hpp>
```

...or include only the pieces you need:

```cpp
#include <vector.hpp>
#include <operations/geometric_product.hpp>
```

Every header is a valid standalone entry point, and that property is maintained
deliberately — each one pulls in its own dependencies. Use the umbrella header
for convenience; use individual headers if compile time matters.

## Your first program

```cpp
#include <iostream>
#include <cliffordcore.hpp>

int main()
{
    CliffordCore::Vector<double> a(1, 2, 3);
    CliffordCore::Vector<double> b(4, 5, 6);

    // The geometric product keeps both the dot and the wedge.
    auto product = a * b;
    std::cout << "dot part:   " << product.scalar.value << "\n";
    std::cout << "wedge part: " << product.bivector.xy << ", "
              << product.bivector.xz << ", " << product.bivector.yz << "\n";

    // A quarter turn about +z, applied to the x axis.
    auto turn = CliffordCore::rotor_from_axis_angle(CliffordCore::Vector<double>(0, 0, 1), 1.5707963);
    auto spun = CliffordCore::rotate(CliffordCore::Vector<double>(1, 0, 0), turn);
    std::cout << "e1 turned:  (" << spun.x << ", " << spun.y << ", " << spun.z << ")\n";
}
```

Compile and run it from the repository root:

```bash
g++ -std=c++17 -Iinclude first.cpp -o first && ./first
```

Expected output:

```text
dot part:   32
wedge part: -3, -6, -3
e1 turned:  (2.67949e-08, 1, 0)
```

`e1` should land exactly on `e2`, so that `2.7e-08` in the x slot is error, not
signal. Most of it comes from the truncated `1.5707963` — a fuller value of
`pi/2` brings it down to around `1e-17`. Rotations leave residue like this
routinely; `examples/print_helpers.hpp` has a `clean()` helper that shows
near-zero components as `0` so the interesting digits stand out.

## Choosing the component type

Every type is a template over an arithmetic `T`, enforced by
`static_assert(std::is_arithmetic<T>::value, ...)`.

- **`double`** — the default choice, and what the tests use.
- **`float`** — works, but rotors drift off unit length noticeably faster. See
  the drift section of `examples/03_rotations.cpp` for measured numbers, and
  renormalise periodically.
- **`long double`** — works.
- **Integer types** — compile, but `norm`, `normalize`, `exp` and `log` all
  truncate. Not useful in practice.

## constexpr

The whole library is marked `constexpr`, but that splits into two groups.

**Portably constant-evaluable** — pure arithmetic, works on any conforming
compiler: the dot, wedge and geometric products, `dual`, `reverse`, the grade
projections, `sandwich`, and mixed-grade `+` and `-`.

```cpp
constexpr CliffordCore::Vector<double> e1(1, 0, 0), e2(0, 1, 0);
static_assert((e1 ^ e2).xy == 1.0, "checked at compile time");
```

**Constant-evaluable only under GCC** — anything reaching `<cmath>`: `norm`,
`normalize`, `exp`, `log`, `rotor_from_axis_angle`, `slerp`, and the
`magnitude()` members. GCC provides constexpr `std::sqrt`, `std::cos` and
friends as an extension; Clang and MSVC do not, so these will not fold at
compile time there. They still work perfectly well at runtime.

## Printing

Every type has `to_string()` and an `operator<<`:

```cpp
std::cout << v << "\n";              // 1*e1 + 2*e2 + 3*e3
std::cout << m << "\n";              // (1) + (2*e1 + ...) + (5*e12 + ...) + (8*e123)
std::string s = r.to_string();
```

Numbers are formatted at round-trip precision, so what you read back parses to
the same value and `float`, `double` and `long double` are distinguishable. That
matters here: this library routinely produces residuals around 1e-16, and a
fixed six-decimal format would print all of them as `0.000000`.

`operator<<` lives in `operations/stream.hpp`, which the umbrella header
includes; the type headers deliberately do not pull in `<ostream>`.

To read components directly instead:

| Type | Fields |
| --- | --- |
| `Scalar<T>` | `.value` |
| `Vector<T>` | `.x` `.y` `.z` |
| `Bivector<T>` | `.xy` `.xz` `.yz` |
| `Trivector<T>` | `.e123` (not `.value`) |
| `Rotor<T>` | `.scalar` and `.bivector` |
| `Multivector<T>` | `.scalar` `.vector` `.bivector` `.trivector` |

`examples/print_helpers.hpp` adds a `clean()` that shows near-zero components as
`0`, which keeps rotation output readable.

## Comparing

`operator==` is **exact**, component by component, and works in a constant
expression. For anything that has been through a rotation or an inverse, use
`approx_equal`:

```cpp
if (approx_equal(rotate(v, r), expected)) { ... }        // default tolerance
if (approx_equal(a, b, 1e-9)) { ... }                    // or state your own
```

The two are separate on purpose: a tolerant `==` is not transitive, which breaks
sorting and associative containers.

## Where next

- **`conventions.md`** — the sign rules. Read this before debugging a sign.
- **`cheatsheet.md`** — return-type matrices and a task-to-call index.
- **`rotations.md`** — rotors in practice.
- **`examples/`** — five runnable programs; `./build_examples.sh` builds and
  runs them all.
