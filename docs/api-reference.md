# API reference

Every entity lives in `namespace CliffordCore`. Every type is
`template<typename T>` with `static_assert(std::is_arithmetic<T>::value, ...)`,
and every function is `constexpr`. `T` is omitted below for readability.

This page is a hand-written catalogue meant to be readable straight from the
repository. For a generated reference with cross-links, run
`cd docs && doxygen Doxyfile`.

---

## Core types

### `cliffordcore.hpp`

Umbrella header. Includes every type and operation; declares nothing itself.

### `scalar.hpp` — `Scalar`

Grade 0. Member: `value`.

| Signature | Returns |
| --- | --- |
| `Scalar()`, `Scalar(T val)` | — |
| `operator+`, `operator-`, `operator*`, `operator/` `(Scalar)` | `Scalar` |
| `operator-()` | `Scalar` |
| `operator*(Vector \| Bivector \| Trivector \| Multivector \| Rotor)` | that same type |
| `operator+=`, `-=`, `*=`, `/=` `(Scalar)` | `Scalar&` |

### `vector.hpp` — `Vector`

Grade 1. Members: `x`, `y`, `z`.

| Signature | Returns |
| --- | --- |
| `Vector()`, `Vector(T x, T y, T z)` | — |
| `magnitude()` | `Scalar` |
| `operator+`, `operator-` `(Vector)`, `operator-()` | `Vector` |
| `operator*`, `operator/` `(Scalar)` | `Vector` |
| `operator+=`, `-=` `(Vector)`; `*=`, `/=` `(T or Scalar)` | `Vector&` |
| free `operator*(T, Vector)` | `Vector` |

### `bivector.hpp` — `Bivector`

Grade 2. Members: `xy`, `xz`, `yz` = `e1e2`, `e1e3`, `e2e3`. Same shape as
`Vector`: `magnitude()`, `+`, `-`, unary `-`, `*`/`/` with `Scalar`, the
compound forms, and free `operator*(T, Bivector)`.

### `trivector.hpp` — `Trivector`

Grade 3, the pseudoscalar. Member: **`e123`** (not `value`).

| Signature | Returns | Note |
| --- | --- | --- |
| `Trivector()`, `Trivector(const T&)` | — | implicit converting ctor |
| `magnitude()` | `Scalar` | |
| `operator+`, `operator-` `(Trivector)`, `operator-()` | `Trivector` | |
| `operator*(Trivector)` | **`Scalar`** | `-ab`, since `e123^2 = -1` |
| `operator/(Trivector)` | **`Scalar`** | `a/b`, the minus signs cancel |
| `operator*`, `operator/` `(T)` | `Trivector` | exact-match, resolves ambiguity |
| `operator*`, `operator/` `(Scalar)` | `Trivector` | scaling |
| `operator+=`, `-=`, `*=`, `/=` | `Trivector&` | |
| free `operator*(T, Trivector)` | `Trivector` | |

### `multivector.hpp` — `Multivector`

All grades. Members: `scalar`, `vector`, `bivector`, `trivector`.

| Signature | Returns |
| --- | --- |
| `Multivector()`, `Multivector(Scalar, Vector, Bivector, Trivector)` | — |
| `Multivector(const Rotor&)` | — (implicit widening) |
| `operator+`, `operator-` `(Multivector)`, `operator-()` | `Multivector` |
| `operator*`, `operator/` `(Scalar)` | `Multivector` |
| `operator+=`, `-=`, `*=`, `/=` | `Multivector&` |
| free `operator*(T, Multivector)` | `Multivector` |

### `rotor.hpp` — `Rotor`

Grades 0 + 2. Members: `scalar`, `bivector`. No `magnitude()` — use `norm(r)`.

| Signature | Returns |
| --- | --- |
| `Rotor()`, `Rotor(Scalar, Bivector)` | — |
| `explicit Rotor(const Multivector&)` | — (**explicit**; drops grades 1 and 3) |
| `operator+`, `operator-` `(Rotor)` | `Rotor` |
| `operator*`, `operator/` `(Scalar)` | `Rotor` |
| `operator+=`, `-=`, `*=`, `/=` | `Rotor&` |
| free `operator*(T, Rotor)` | `Rotor` |

---

## Products

### `dot_product.hpp`

- `operator|(Vector, Vector)` → `Scalar`
- `dot_product(Vector, Vector)` → `Scalar`

Vector-vector only; `contraction.hpp` extends `operator|` across grades.

### `wedge_product.hpp`

- `operator^(Vector, Vector)` → `Bivector`
- `wedge_product(Vector, Vector)` → `Bivector`

Grade 1 ∧ grade 2 lives in `mixed_products.hpp`.

### `geometric_product.hpp`

- `geometric_product(Vector, Vector)` → `Multivector`
- `operator*(Vector, Vector)` → `Multivector`
- `rotor_product(Vector, Vector)` → `Rotor`
- `geometric_product(Multivector, Multivector)` → `Multivector` — **the full
  8x8 multiplication table; every mixed product routes through this**
- `operator*(Multivector, Multivector)` → `Multivector`
- `operator*(Rotor, Rotor)` → `Rotor`
- `rotor_product(Rotor, Rotor)` → `Rotor`

### `mixed_products.hpp`

`operator*` for every remaining pair, all returning `Multivector`: vector with
bivector, bivector with bivector, vector with trivector, bivector with trivector,
rotor with vector, and multivector with each of vector, bivector, trivector and
rotor — in both operand orders.

Also the wedge beyond grade 1:

- `operator^(Vector, Bivector)` → `Trivector`
- `operator^(Bivector, Vector)` → `Trivector`
- `wedge_product(...)` for both

With these, `operator*` is defined for **every** operand pair; see
`cheatsheet.md` for the full matrix.

---

## Sums

### `addition.hpp`

Twenty `operator+` overloads covering every mixed-grade pair in both orders, all
returning `Multivector`. Plus:

- `rotor_sum(Scalar, Bivector)` → `Rotor`
- `to_rotor(Multivector)` → `Rotor` (drops grades 1 and 3)
- `to_multivector(Rotor)` → `Multivector`

Internal: `detail::promote(X)` → `Multivector`, for each grade.

### `subtraction.hpp`

The same twenty pairs as `operator-`. Reuses `detail::promote`.

---

## Magnitudes and involutions

### `norm.hpp`

- `norm(x)` → `Scalar`, for `Scalar`, `Vector`, `Bivector`, `Trivector`,
  `Multivector`, `Rotor`
- `squared_norm(x)` → `Scalar`, for all six types

### `normalize.hpp`

- `normalize(x)` → same type, for all six types. `normalize(Scalar)` gives the
  sign, +1 or -1.

A zero input is returned unchanged rather than producing `NaN`.

### `reverse.hpp`

- `reverse(x)` → same type, for all six. Grades 0 and 1 keep their sign,
  grades 2 and 3 are negated.

### `inverse.hpp`

- `inverse(x)` → same type, for `Scalar`, `Vector`, `Bivector`, `Trivector`,
  `Multivector`, `Rotor`
`inverse(Scalar)` does not guard against zero. `inverse(Multivector)` returns a
zero multivector when the input is not invertible. See `conventions.md` for the
sign of each.

---

## Structure

### `dual.hpp`

- `dual(Scalar)` → `Trivector`
- `dual(Vector)` → `Bivector`
- `dual(Bivector)` → `Vector`
- `dual(Trivector)` → `Scalar`

Not defined for `Multivector` or `Rotor`.

### `grade.hpp`

- `grade0(Multivector)` → `Scalar`
- `grade1(Multivector)` → `Vector`
- `grade2(Multivector)` → `Bivector`
- `grade3(Multivector)` → `Trivector`
- `grade0(Rotor)` → `Scalar`
- `grade2(Rotor)` → `Bivector`

Named per grade rather than `grade<N>()` because each returns a different type.

### `involutions.hpp`

- `involute(x)` → same type, all six — grade `k` scaled by `(-1)^k`
- `conjugate(x)` → same type, all six — grade `k` scaled by `(-1)^(k(k+1)/2)`

Together with `reverse` these are the three involutions; `conjugate` equals
`reverse` composed with `involute`.

### `contraction.hpp`

- `left_contraction(a, b)`, also spelled `a | b` — for the blade pairs where the
  result grade is non-negative: `V,V`→`Scalar`; `V,B`→`Vector`;
  `V,T`→`Bivector`; `B,B`→`Scalar`; `B,T`→`Vector`; `T,T`→`Scalar`
- `right_contraction(a, b)` — the reverse-mirror, so the operands swap roles
- `scalar_product(a, b)` → `Scalar` — the grade 0 part of any product

`operator|` on two vectors is unchanged; this only extends it across grades.

### `geometry.hpp`

- `reflect(Vector, Vector)` → `Vector` — reflect in the plane normal to the
  second argument
- `reflect(Vector, Bivector)` → `Vector` — reflect in the plane itself
- `project(Vector, Vector | Bivector)` → `Vector`
- `reject(Vector, Vector | Bivector)` → `Vector`

None require a unit second argument; each divides through by `inverse`.
`project(v, x) + reject(v, x) == v`.

### `comparison.hpp`

- `operator==`, `operator!=` — **exact**, component by component, `constexpr`
- `approx_equal(a, b, tolerance)` — tolerance defaults to 100 epsilons

Separate on purpose: a tolerant `==` is not transitive.

### `stream.hpp`

- `operator<<(std::ostream&, x)` for all six, delegating to `to_string()`

Kept apart so the type headers need not include `<ostream>`.

---

## Rotations

### `exp.hpp`

- `exp(Bivector)` → `Rotor` — returns the identity rotor for a zero bivector

### `log.hpp`

- `log(Rotor)` → `Bivector` — clamps the scalar part into `[-1, 1]` before
  `std::acos`; returns a zero bivector when the bivector part is zero

Both collide by name with `<cmath>`, so qualify them as `CliffordCore::exp` /
`CliffordCore::log` when `using namespace` is in play.

### `sandwich.hpp`

- `sandwich(x, Rotor)` → same type as `x`, for `Vector`, `Bivector`,
  `Trivector` and `Multivector` — computes `R x reverse(R)`
- `rotate(x, Rotor)` — alias for each

**Object first, rotor second.** Rotating a `Trivector` is the identity for a
unit rotor: the pseudoscalar is central.

### `rotor_construction.hpp`

- `identity_rotor<T>()` → `Rotor` — `T` must be given explicitly
- `rotor_from_axis_angle(Vector axis, T angle)` → `Rotor`
- `rotor_between(Vector from, Vector to)` → `Rotor`
- `slerp(Rotor from, Rotor to, T t)` → `Rotor`
