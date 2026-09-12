# Cheat sheet

One screen. Everything here is verified by compilation and by
`tests/test_core.cpp`.

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
