# API reference

**This page catalogues Cl(3,0)**, `namespace CliffordCore::Cl3`; the snippets
below assume `namespace ga = CliffordCore::Cl3;`. Every type is
`template<typename T>` with `static_assert(std::is_arithmetic<T>::value, ...)`,
and every function is `constexpr`. `T` is omitted below for readability.

Cl(3,0,1) mirrors this structure one directory over, in
`cliffordcore/pga/` and `namespace CliffordCore::PGA`, with the same header
names and the same operation names. What differs is the *meaning* of each type
and several of the formulas; [pga.md](pga.md) documents those, and the PGA
section at the end of this page lists what the namespace adds and where the
signatures differ.

This page is a hand-written catalogue meant to be readable straight from the
repository. For a generated reference with cross-links, run
`cd docs && doxygen Doxyfile`.

---

## Core types

### `cliffordcore/cl3.hpp`

Umbrella for Cl(3,0). Includes every type and operation below; declares nothing
itself. Headers are listed by their name within `cliffordcore/cl3/`.

### `cliffordcore.hpp`

**Not an umbrella.** Holds only `CLIFFORDCORE_VERSION_MAJOR` / `_MINOR` /
`_PATCH` / `_STRING` and the Doxygen mainpage. There is deliberately no header
that pulls in every algebra.

### `cliffordcore/detail/`

Algebra-independent helpers in `CliffordCore::detail`, shared by every algebra:
`format_component` (round-trip number formatting, in `format.hpp`) and
`default_tolerance` / `close` (in `compare.hpp`). Not part of the public API.

Note that each algebra also has its own `detail` — `CliffordCore::Cl3::detail`
holds `promote`, which widens to *that algebra's* `Multivector`. The inner one
shadows the outer, so library code spells the shared helpers in full.

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

**A zero input is a precondition violation, not a value.** Every overload
asserts that the norm is non-zero, since a zero element has no direction to
preserve. The assert compiles out under `NDEBUG`, where the argument is returned
unchanged instead — a release build must not divide by zero and hand back `NaN`.

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

Both collide by name with `<cmath>`, so qualify them as `ga::exp` / `ga::log`
when `using namespace` is in play. The alias makes that the natural spelling
anyway.

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

---

## Cl(3,0,1) — `namespace CliffordCore::PGA`

Header `cliffordcore/pga.hpp`, directory `cliffordcore/pga/`. Same layout, same
operation names, same `constexpr`-everywhere rule. Read [pga.md](pga.md) for the
conventions; this section is the delta from the Cl(3,0) catalogue above.

### Types

| Type | Grade | Members | Geometry |
| --- | --- | --- | --- |
| `Scalar` | 0 | `value` | |
| `Vector`, `Plane` | 1 | `e0 e1 e2 e3` | a plane |
| `Bivector`, `Line`, `Twist` | 2 | `e01 e02 e03 e12 e13 e23` | a line |
| `Trivector`, `Point` | 3 | `e012 e013 e023 e123` | a point |
| `Quadvector` | 4 | `e0123` | the pseudoscalar |
| `Multivector` | mixed | the five parts above | |
| `Rotor` | 0+2 | `scalar, e12, e13, e23` | rotation about an origin axis |
| `Translator` | 0+2 | `scalar, e01, e02, e03` | translation |
| `Motor` | 0+2+4 | `scalar, bivector, quadvector` | any rigid motion |

`Rotor`, `Translator` and `Motor` widen implicitly to `Multivector`, and the
first two to `Motor`. Every narrowing constructor is `explicit`; the named
forms are `to_rotor`, `to_translator`, `to_motor`, `to_multivector`.

### What the shared headers do differently

- `geometric_product.hpp` — the 16x16 table. Adds `motor_product(Vector, Vector)`
  and `translator_product`; `Rotor * Rotor`, `Translator * Translator` and
  `Motor * Motor` all keep their own type.
- `mixed_products.hpp` — every pair of even types returns `Motor`, everything
  else `Multivector`. **`Trivector * Trivector` is here and returns a
  `Multivector`** (grades 0 and 2), not a `Scalar` as in Cl(3,0).
- `norm.hpp` — adds `ideal_norm(x)` for every type. `norm` sees only the
  components without `e0`, so every ideal element has zero norm.
- `normalize.hpp` — no `Quadvector` overload. `normalize(Motor)` divides by the
  square root of a study number, not by the norm.
- `inverse.hpp` — no `Quadvector` overload. `Bivector` and `Motor` carry a study
  number correction; `Multivector` uses the Hitzer–Sangwine n = 4 formula.
- `dual.hpp` — `dual` is the right complement, not a product with the
  pseudoscalar, and `undual` is its exact inverse. Both cover all five grades
  and `Multivector`.
- `grade.hpp` — adds `grade4`, and `grade0`/`grade2` for `Translator`, plus
  `grade0`/`grade2`/`grade4` for `Motor`.
- `contraction.hpp` — the ten non-decreasing pairs among the five blade types.
- `exp.hpp` / `log.hpp` — `exp(Bivector)` → `Motor`;
  `log(Motor | Rotor | Translator)` → `Bivector`.
- `sandwich.hpp` — `sandwich(x, Motor | Rotor | Translator)` for `Vector`,
  `Bivector`, `Trivector` and `Multivector`, plus `rotate`, `translate` and
  `transform` as the named forms.

### Headers with no Cl(3,0) counterpart

`regressive_product.hpp`

- `regressive_product(a, b)`, `join(a, b)` — for `(Point, Point)` → `Line`,
  `(Point, Line)` and `(Line, Point)` → `Plane`, `(Line, Line)` → `Scalar`
- `meet(a, b)` — the wedge under its geometric name

`primitives.hpp`

- `plane(a, b, c, d)`, `point(x, y, z)`, `ideal_point(dx, dy, dz)`,
  `origin<T>()`, `ideal_plane<T>()`
- `line_through_origin(dx, dy, dz)`, `line_through_points(P, Q)`,
  `line_from_planes(p, q)`, `plane_through_points(P, Q, R)`
- `position(Point)`, `direction(Point)`, `direction(Line)`, `moment(Line)`,
  `normal(Plane)`, `offset(Plane)` — all `std::array<T, 3>` except `offset`
- `is_ideal(x)` for `Plane`, `Line` and `Point`

`motor_construction.hpp` (`rotor_construction.hpp`'s counterpart)

- `identity_rotor<T>()`, `identity_translator<T>()`, `identity_motor<T>()`
- `translator(dx, dy, dz)`, `translator_between(Point, Point)`
- `rotor_from_axis_angle(dx, dy, dz, angle)` — note the **direction as three
  numbers**, not a `Vector`, since a PGA `Vector` is a plane
- `motor_from_line_angle(Line, angle)`, `screw(Line, angle, distance)`
- `motor_between(Plane, Plane)`, `motor_from_rotor_translator(Rotor, Translator)`
- `slerp` for `Motor`, `Rotor` and `Translator`

`geometry.hpp` (same name, different contents)

- `reflect(x, Plane)` for `Plane`, `Line` and `Point`
- `project(Point, Plane)`, `project(Point, Line)`, `project(Line, Plane)`
- `distance(Point, Point)`, `distance(Point, Plane)` (signed),
  `distance(Point, Line)`
- `angle(Plane, Plane)`, `angle(Line, Line)`

`distance` and `angle` return a raw `T` and are `inline`, not `constexpr` —
they use `std::acos` and `std::sqrt`. There is no `reject`.
