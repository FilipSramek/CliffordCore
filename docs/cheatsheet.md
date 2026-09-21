# Cheat sheet

One screen, for Cl(3,0). Everything here is verified by compilation and by
`tests/test_core.cpp`. The Cl(3,0,1) equivalent is at the bottom.

```cpp
#include <cliffordcore/cl3.hpp>
namespace ga = CliffordCore::Cl3;
```

## `operator*` return types

| `*` | Scalar | Vector | Bivector | Trivector | Multivector | Rotor |
| --- | --- | --- | --- | --- | --- | --- |
| **Scalar** | Scalar | Vector | Bivector | Trivector | Multivector | Rotor |
| **Vector** | Vector | Multivector | Multivector | Multivector | Multivector | Multivector |
| **Bivector** | Bivector | Multivector | Multivector | Multivector | Multivector | Multivector |
| **Trivector** | Trivector | Multivector | Multivector | **Scalar** | Multivector | Multivector |
| **Multivector** | Multivector | Multivector | Multivector | Multivector | Multivector | Multivector |
| **Rotor** | Rotor | Multivector | Multivector | Multivector | Multivector | Rotor |

Every pair is defined. Two cells to remember: `Trivector * Trivector` is a **Scalar** (`e123^2 = -1`),
and `Rotor * Rotor` stays a **Rotor** so composition never widens.

## `operator+` and `operator-` return types

Every mixed pair gives `Multivector`. Same-grade pairs keep their type.
`operator-` has exactly the same shape.

| `+` `-` | Scalar | Vector | Bivector | Trivector | Multivector | Rotor |
| --- | --- | --- | --- | --- | --- | --- |
| **Scalar** | Scalar | Multivector | Multivector | Multivector | Multivector | — |
| **Vector** | Multivector | Vector | Multivector | Multivector | Multivector | — |
| **Bivector** | Multivector | Multivector | Bivector | Multivector | Multivector | — |
| **Trivector** | Multivector | Multivector | Multivector | Trivector | Multivector | — |
| **Multivector** | Multivector | Multivector | Multivector | Multivector | Multivector | Multivector |
| **Rotor** | — | — | — | — | — | Rotor |

### The rotor row

`Rotor` adds only with another `Rotor`, or with a `Multivector` on the left.
For anything else, widen it first:

```cpp
auto sum = to_multivector(r) + v;   // Multivector
```

This used to be worse. `Rotor` once had an *implicit* constructor from
`Multivector`, so `r + m` compiled, narrowed `m` down to a rotor, and threw
away grades 1 and 3 with no diagnostic — while `m + r` kept everything. The
constructor is `explicit` now, so the lossy spelling is a compile error rather
than a silent wrong answer.

## Common tasks

| To do this | Call this |
| --- | --- |
| dot product / left contraction | `a \| b` or `left_contraction(a, b)` |
| right contraction | `right_contraction(a, b)` |
| grade 0 part of a product | `scalar_product(a, b)` |
| wedge / outer product | `a ^ b` or `wedge_product(a, b)` |
| geometric product | `a * b` or `geometric_product(a, b)` |
| get a rotor from two vectors | `rotor_product(a, b)` |
| length | `norm(x)` |
| length squared (no sqrt) | `squared_norm(x)` |
| unit-length version | `normalize(x)` |
| flip grades 2 and 3 | `reverse(x)` |
| flip grades 1 and 2 | `conjugate(x)` |
| flip the odd grades, 1 and 3 | `involute(x)` |
| multiplicative inverse | `inverse(x)` |
| swap grade k with grade 3-k | `dual(x)` |
| pull one grade out | `grade0(m)` … `grade3(m)` |
| bivector to rotor | `exp(b)` |
| rotor to bivector | `log(r)` |
| no-op rotor | `identity_rotor<double>()` |
| rotor from axis and angle | `rotor_from_axis_angle(axis, angle)` |
| rotor taking one vector to another | `rotor_between(from, to)` |
| apply a rotation | `rotate(v, r)` or `sandwich(v, r)` |
| rotate a plane or a multivector | `rotate(b, r)`, `rotate(m, r)` |
| compose two rotations | `r1 * r2` |
| undo a rotation | `inverse(r)`, or `reverse(r)` if unit |
| interpolate two rotations | `slerp(from, to, t)` |
| angle between vectors | `atan2(norm(a ^ b).value, (a \| b).value)` |
| signed volume | `((a ^ b) ^ c).e123` |
| triangle area | `norm(a ^ b).value / 2` |
| reflect v in the plane normal to n | `reflect(v, n)` |
| reflect v in the plane B | `reflect(v, B)` |
| project v onto a line, or a plane | `project(v, u)`, `project(v, B)` |
| the leftover perpendicular part | `reject(v, u)`, `reject(v, B)` |
| compare exactly | `a == b`, `a != b` |
| compare within a tolerance | `approx_equal(a, b, tol)` |
| print | `os << a`, or `a.to_string()` |

## Printing

`operator<<` and `to_string()` work on every type, at round-trip precision:

```cpp
std::cout << v << std::endl;     // 1*e1 + 2*e2 + 3*e3
std::string s = m.to_string();
```

To read components directly instead:

| Type | Fields |
| --- | --- |
| `Scalar<T>` | `.value` |
| `Vector<T>` | `.x` `.y` `.z` |
| `Bivector<T>` | `.xy` `.xz` `.yz` |
| `Trivector<T>` | `.e123` |
| `Rotor<T>` | `.scalar.value`, `.bivector.xy` … |
| `Multivector<T>` | `.scalar` `.vector` `.bivector` `.trivector` |

## Signs, in four lines

```text
e1^2 = e2^2 = e3^2 = +1        e12^2 = e13^2 = e23^2 = e123^2 = -1
bivector basis is (xy, xz, yz) = (e1e2, e1e3, e2e3)     -- xz, not zx
dual(dual(a)) == -a            exp(theta*B) rotates by 2*theta
sandwich(v, r) = R v ~R        e12 rotates about z, not about x
```

---

## Cl(3,0,1), on one screen

```cpp
#include <cliffordcore/pga.hpp>
namespace ga = CliffordCore::PGA;
```

Verified by `tests/test_pga.cpp`. Full conventions in [pga.md](pga.md).

### What each grade is

| Type | Is | Alias |
| --- | --- | --- |
| `Vector` | a **plane** `a x + b y + c z + d = 0` | `Plane` |
| `Bivector` | a **line** | `Line`, `Twist` |
| `Trivector` | a **point** | `Point` |
| `Quadvector` | the pseudoscalar `e0123` | |
| `Rotor` / `Translator` / `Motor` | rotation about an origin axis / translation / any rigid motion | |

### Building things

```cpp
auto p  = ga::plane(a, b, c, d);              // a x + b y + c z + d = 0
auto P  = ga::point(x, y, z);                 // weight 1
auto O  = ga::origin<double>();
auto L  = ga::line_through_points(P, Q);      // or line_from_planes(p, q)
auto ax = ga::line_through_origin(dx, dy, dz);
```

### Moving things

```cpp
auto t = ga::translator(dx, dy, dz);
auto r = ga::rotor_from_axis_angle(dx, dy, dz, angle);   // axis through the origin
auto m = ga::motor_from_line_angle(L, angle);            // ANY line
auto s = ga::screw(L, angle, distance);                  // rotate and slide

auto moved = ga::sandwich(P, m);     // works for planes, lines, points, multivectors
auto back  = ga::sandwich(moved, ga::reverse(m));
auto half  = ga::slerp(identity, m, 0.5);
```

### Asking questions

```cpp
ga::meet(p, q)            // plane ^ plane -> the line they share
ga::join(P, Q)            // point v point -> the line through both
(p ^ P).e0123             // zero iff the point lies on the plane
ga::distance(P, Q)        // also (P, plane) -- signed -- and (P, line)
ga::angle(p, q)           // also (line, line)
ga::project(P, p)         // also (P, line) and (line, plane)
ga::reflect(P, mirror)
ga::position(P)           // std::array<T,3>; direction(), normal(), moment()
ga::is_ideal(x)           // at infinity?
```

### Components

| Type | Fields |
| --- | --- |
| `Vector` | `.e0` `.e1` `.e2` `.e3` |
| `Bivector` | `.e01` `.e02` `.e03` `.e12` `.e13` `.e23` |
| `Trivector` | `.e012` `.e013` `.e023` `.e123` |
| `Quadvector` | `.e0123` |
| `Rotor` | `.scalar.value` `.e12` `.e13` `.e23` |
| `Translator` | `.scalar.value` `.e01` `.e02` `.e03` |
| `Motor` | `.scalar` `.bivector` `.quadvector` |
| `Multivector` | `.scalar` `.vector` `.bivector` `.trivector` `.quadvector` |

### Signs and traps, in six lines

```text
e0^2 = 0, e1^2 = e2^2 = e3^2 = +1        e123^2 = -1,  e0123^2 = 0
point(x,y,z) = e123 - x e023 + y e013 - z e012      x,y,z axes = e23, -e13, e12
dual is the COMPLEMENT, not *e0123       dual(dual(a)) == -a on grades 1 and 3
Trivector * Trivector is a Multivector   (grades 0 and 2) -- NOT a Scalar as in Cl3
norm() ignores every e0 component        ideal_norm() is the other half
a Motor's inverse/normalize need the study number; reverse(m)/|m|^2 is not enough
```
