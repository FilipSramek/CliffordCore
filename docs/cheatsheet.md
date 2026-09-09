# Cheat sheet

One screen. Everything here is verified by compilation and by
`tests/test_core.cpp`.

## `operator*` return types

| `*` | Scalar | Vector3 | Bivector3 | Trivector3 | Multivector3 | Rotor3 |
| --- | --- | --- | --- | --- | --- | --- |
| **Scalar** | Scalar | Vector3 | Bivector3 | Trivector3 | Multivector3 | Rotor3 |
| **Vector3** | Vector3 | Multivector3 | Multivector3 | Multivector3 | Multivector3 | Multivector3 |
| **Bivector3** | Bivector3 | Multivector3 | Multivector3 | Multivector3 | Multivector3 | — |
| **Trivector3** | Trivector3 | Multivector3 | Multivector3 | **Scalar** | Multivector3 | — |
| **Multivector3** | Multivector3 | Multivector3 | Multivector3 | Multivector3 | Multivector3 | Multivector3 |
| **Rotor3** | Rotor3 | Multivector3 | — | — | Multivector3 | Rotor3 |

`—` marks a combination that does not compile. Work around it with
`to_multivector(r)`.

Two cells to remember: `Trivector3 * Trivector3` is a **Scalar** (`e123^2 = -1`),
and `Rotor3 * Rotor3` stays a **Rotor3** so composition never widens.

## `operator+` and `operator-` return types

Every mixed pair gives `Multivector3`. Same-grade pairs keep their type.
`operator-` has exactly the same shape.

| `+` `-` | Scalar | Vector3 | Bivector3 | Trivector3 | Multivector3 | Rotor3 |
| --- | --- | --- | --- | --- | --- | --- |
| **Scalar** | Scalar | Multivector3 | Multivector3 | Multivector3 | Multivector3 | — |
| **Vector3** | Multivector3 | Vector3 | Multivector3 | Multivector3 | Multivector3 | — |
| **Bivector3** | Multivector3 | Multivector3 | Bivector3 | Multivector3 | Multivector3 | — |
| **Trivector3** | Multivector3 | Multivector3 | Multivector3 | Trivector3 | Multivector3 | — |
| **Multivector3** | Multivector3 | Multivector3 | Multivector3 | Multivector3 | Multivector3 | Multivector3 |
| **Rotor3** | — | — | — | — | **Rotor3 (lossy!)** | Rotor3 |

### Watch the rotor row

`Rotor3 + Multivector3` and `Multivector3 + Rotor3` are **not the same
operation**, and one of them silently loses data:

```cpp
Rotor3<double>       r(Scalar<double>(1), Bivector3<double>(1, 0, 0));
Multivector3<double> m(Scalar<double>(10), Vector3<double>(20, 30, 40),
                       Bivector3<double>(50, 0, 0), Trivector3<double>(60));

auto a = r + m;   // Rotor3       -- vector (20,30,40) and trivector 60 DISCARDED
auto b = m + r;   // Multivector3 -- everything preserved
```

`Rotor3` has an implicit constructor from `Multivector3`, so `r + m` converts
`m` down to a rotor first and throws away grades 1 and 3 without a warning.
**Put the multivector on the left**, or convert explicitly with
`to_multivector(r)`.

The other rotor cells are genuinely undefined: use `to_multivector(r)`.

## Common tasks

| To do this | Call this |
| --- | --- |
| dot product | `a \| b` or `dot_product(a, b)` |
| wedge / outer product | `a ^ b` or `wedge_product(a, b)` |
| geometric product | `a * b` or `geometric_product(a, b)` |
| get a rotor from two vectors | `rotor_product(a, b)` |
| length | `norm(x)` |
| length squared (no sqrt) | `squared_norm(x)` |
| unit-length version | `normalize(x)` |
| flip grades 2 and 3 | `reverse(x)` |
| flip grades 1 and 2 | `conjugate(m)` |
| multiplicative inverse | `inverse(x)` |
| swap grade k with grade 3-k | `dual(x)` |
| pull one grade out | `grade0(m)` … `grade3(m)` |
| bivector to rotor | `exp(b)` |
| rotor to bivector | `log(r)` |
| no-op rotor | `identity_rotor<double>()` |
| rotor from axis and angle | `rotor_from_axis_angle(axis, angle)` |
| rotor taking one vector to another | `rotor_between(from, to)` |
| apply a rotation | `rotate(v, r)` or `sandwich(v, r)` |
| compose two rotations | `r1 * r2` |
| undo a rotation | `inverse(r)`, or `reverse(r)` if unit |
| interpolate two rotations | `slerp(from, to, t)` |
| angle between vectors | `atan2(norm(a ^ b).value, (a \| b).value)` |
| signed volume | `((a ^ b) ^ c).e123` |
| triangle area | `norm(a ^ b).value / 2` |
| project v onto unit u | `u * (v \| u)` |
| reflect v in plane with unit normal n | `grade1(-(n * v * n))` |

## Reading components

There is no `operator<<`, so print these directly:

| Type | Fields |
| --- | --- |
| `Scalar<T>` | `.value` |
| `Vector3<T>` | `.x` `.y` `.z` |
| `Bivector3<T>` | `.xy` `.xz` `.yz` |
| `Trivector3<T>` | `.e123` |
| `Rotor3<T>` | `.scalar.value`, `.bivector.xy` … |
| `Multivector3<T>` | `.scalar` `.vector` `.bivector` `.trivector` |

## Signs, in four lines

```text
e1^2 = e2^2 = e3^2 = +1        e12^2 = e13^2 = e23^2 = e123^2 = -1
bivector basis is (xy, xz, yz) = (e1e2, e1e3, e2e3)     -- xz, not zx
dual(dual(a)) == -a            exp(theta*B) rotates by 2*theta
sandwich(v, r) = R v ~R        e12 rotates about z, not about x
```
