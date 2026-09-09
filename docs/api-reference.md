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
| `operator*(Vector3 \| Bivector3 \| Trivector3 \| Multivector3 \| Rotor3)` | that same type |
| `operator+=`, `-=`, `*=`, `/=` `(Scalar)` | `Scalar&` |

### `vector3.hpp` — `Vector3`

Grade 1. Members: `x`, `y`, `z`.

| Signature | Returns |
| --- | --- |
| `Vector3()`, `Vector3(T x, T y, T z)` | — |
| `magnitude()` | `Scalar` |
| `operator+`, `operator-` `(Vector3)`, `operator-()` | `Vector3` |
| `operator*`, `operator/` `(Scalar)` | `Vector3` |
| `operator+=`, `-=` `(Vector3)`; `*=`, `/=` `(T or Scalar)` | `Vector3&` |
| free `operator*(T, Vector3)` | `Vector3` |

### `bivector3.hpp` — `Bivector3`

Grade 2. Members: `xy`, `xz`, `yz` = `e1e2`, `e1e3`, `e2e3`. Same shape as
`Vector3`: `magnitude()`, `+`, `-`, unary `-`, `*`/`/` with `Scalar`, the
compound forms, and free `operator*(T, Bivector3)`.

### `trivector3.hpp` — `Trivector3`

Grade 3, the pseudoscalar. Member: **`e123`** (not `value`).

| Signature | Returns | Note |
| --- | --- | --- |
| `Trivector3()`, `Trivector3(const T&)` | — | implicit converting ctor |
| `magnitude()` | `Scalar` | |
| `operator+`, `operator-` `(Trivector3)`, `operator-()` | `Trivector3` | |
| `operator*(Trivector3)` | **`Scalar`** | `-ab`, since `e123^2 = -1` |
| `operator/(Trivector3)` | **`Scalar`** | `a/b`, the minus signs cancel |
| `operator*`, `operator/` `(T)` | `Trivector3` | exact-match, resolves ambiguity |
| `operator*`, `operator/` `(Scalar)` | `Trivector3` | scaling |
| `operator+=`, `-=`, `*=`, `/=` | `Trivector3&` | |
| free `operator*(T, Trivector3)` | `Trivector3` | |

### `multivector3.hpp` — `Multivector3`

All grades. Members: `scalar`, `vector`, `bivector`, `trivector`.

| Signature | Returns |
| --- | --- |
| `Multivector3()`, `Multivector3(Scalar, Vector3, Bivector3, Trivector3)` | — |
| `Multivector3(const Rotor3&)` | — (implicit widening) |
| `operator+`, `operator-` `(Multivector3)`, `operator-()` | `Multivector3` |
| `operator*`, `operator/` `(Scalar)` | `Multivector3` |
| `operator+=`, `-=`, `*=`, `/=` | `Multivector3&` |
| free `operator*(T, Multivector3)` | `Multivector3` |

### `rotor3.hpp` — `Rotor3`

Grades 0 + 2. Members: `scalar`, `bivector`. No `magnitude()` — use `norm(r)`.

| Signature | Returns |
| --- | --- |
| `Rotor3()`, `Rotor3(Scalar, Bivector3)` | — |
| `Rotor3(const Multivector3&)` | — (implicit, **drops grades 1 and 3**) |
| `operator+`, `operator-` `(Rotor3)` | `Rotor3` |
| `operator*`, `operator/` `(Scalar)` | `Rotor3` |
| `operator+=`, `-=`, `*=`, `/=` | `Rotor3&` |
| free `operator*(T, Rotor3)` | `Rotor3` |

---

## Products

### `dot_product.hpp`

- `operator|(Vector3, Vector3)` → `Scalar`
- `dot_product(Vector3, Vector3)` → `Scalar`

Vector-vector only. There are no contractions for other grades.

### `wedge_product.hpp`

- `operator^(Vector3, Vector3)` → `Bivector3`
- `wedge_product(Vector3, Vector3)` → `Bivector3`

Grade 1 ∧ grade 2 lives in `mixed_products.hpp`.

### `geometric_product.hpp`

- `geometric_product(Vector3, Vector3)` → `Multivector3`
- `operator*(Vector3, Vector3)` → `Multivector3`
- `rotor_product(Vector3, Vector3)` → `Rotor3`
- `geometric_product(Multivector3, Multivector3)` → `Multivector3` — **the full
  8x8 multiplication table; every mixed product routes through this**
- `operator*(Multivector3, Multivector3)` → `Multivector3`
- `operator*(Rotor3, Rotor3)` → `Rotor3`
- `rotor_product(Rotor3, Rotor3)` → `Rotor3`

### `mixed_products.hpp`

`operator*` for every remaining pair, all returning `Multivector3`: vector with
bivector, bivector with bivector, vector with trivector, bivector with trivector,
rotor with vector, and multivector with each of vector, bivector, trivector and
rotor — in both operand orders.

Also the wedge beyond grade 1:

- `operator^(Vector3, Bivector3)` → `Trivector3`
- `operator^(Bivector3, Vector3)` → `Trivector3`
- `wedge_product(...)` for both

See `cheatsheet.md` for the four undefined cells.

---

## Sums

### `addition.hpp`

Twenty `operator+` overloads covering every mixed-grade pair in both orders, all
returning `Multivector3`. Plus:

- `rotor_sum(Scalar, Bivector3)` → `Rotor3`
- `to_rotor(Multivector3)` → `Rotor3` (drops grades 1 and 3)
- `to_multivector(Rotor3)` → `Multivector3`

Internal: `detail::promote(X)` → `Multivector3`, for each grade.

### `subtraction.hpp`

The same twenty pairs as `operator-`. Reuses `detail::promote`.

---

## Magnitudes and involutions

### `norm.hpp`

- `norm(x)` → `Scalar`, for `Scalar`, `Vector3`, `Bivector3`, `Trivector3`,
  `Multivector3`, `Rotor3`
- `squared_norm(x)` → `Scalar`, for all of those **except `Scalar`**

### `normalize.hpp`

- `normalize(x)` → same type, for `Vector3`, `Bivector3`, `Trivector3`,
  `Multivector3`, `Rotor3`. Not defined for `Scalar`.

A zero input is returned unchanged rather than producing `NaN`.

### `reverse.hpp`

- `reverse(x)` → same type, for `Vector3` (identity), `Bivector3` (negated),
  `Trivector3` (negated), `Multivector3`, `Rotor3`. Not defined for `Scalar`.

### `inverse.hpp`

- `inverse(x)` → same type, for `Scalar`, `Vector3`, `Bivector3`, `Trivector3`,
  `Multivector3`, `Rotor3`
- `conjugate(Multivector3)` → `Multivector3` — Clifford conjugation

`inverse(Scalar)` does not guard against zero. `inverse(Multivector3)` returns a
zero multivector when the input is not invertible. See `conventions.md` for the
sign of each.

---

## Structure

### `dual.hpp`

- `dual(Scalar)` → `Trivector3`
- `dual(Vector3)` → `Bivector3`
- `dual(Bivector3)` → `Vector3`
- `dual(Trivector3)` → `Scalar`

Not defined for `Multivector3` or `Rotor3`.

### `grade.hpp`

- `grade0(Multivector3)` → `Scalar`
- `grade1(Multivector3)` → `Vector3`
- `grade2(Multivector3)` → `Bivector3`
- `grade3(Multivector3)` → `Trivector3`
- `grade0(Rotor3)` → `Scalar`
- `grade2(Rotor3)` → `Bivector3`

Named per grade rather than `grade<N>()` because each returns a different type.

---

## Rotations

### `exp.hpp`

- `exp(Bivector3)` → `Rotor3` — returns the identity rotor for a zero bivector

### `log.hpp`

- `log(Rotor3)` → `Bivector3` — clamps the scalar part into `[-1, 1]` before
  `std::acos`; returns a zero bivector when the bivector part is zero

Both collide by name with `<cmath>`, so qualify them as `CliffordCore::exp` /
`CliffordCore::log` when `using namespace` is in play.

### `sandwich.hpp`

- `sandwich(Vector3, Rotor3)` → `Vector3` — computes `R v reverse(R)`
- `rotate(Vector3, Rotor3)` → `Vector3` — alias

**Vector first, rotor second.**

### `rotor_construction.hpp`

- `identity_rotor<T>()` → `Rotor3` — `T` must be given explicitly
- `rotor_from_axis_angle(Vector3 axis, T angle)` → `Rotor3`
- `rotor_between(Vector3 from, Vector3 to)` → `Rotor3`
- `slerp(Rotor3 from, Rotor3 to, T t)` → `Rotor3`
