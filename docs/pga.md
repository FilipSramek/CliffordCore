# PGA: Cl(3,0,1)

`CliffordCore::PGA` is the projective geometric algebra of 3D Euclidean space.
Where `Cl3` gives you directions and rotations about the origin, PGA gives you
planes, lines and points **anywhere in space**, plus the rigid motions between
them: translations, rotations about arbitrary axes, and the screws that combine
both.

It reuses the Cl3 type names, so switching algebra is one line:

```cpp
#include <cliffordcore/pga.hpp>
namespace ga = CliffordCore::PGA;
```

**It reuses the names, not the meanings.** A `Vector` here is a plane, not a
direction; a `Trivector` is a point, not the pseudoscalar. None of the sign
conventions in [conventions.md](conventions.md) carry over. Everything below
was verified numerically against an independent model of the algebra before it
was written down.

## The basis

Four basis vectors, `e0 e1 e2 e3`, with

```
e0^2 = 0        e1^2 = e2^2 = e3^2 = +1
```

`e0` is the degenerate one. It is what carries the *offset* of an object from
the origin, and its square being zero is why translations work out as versors
at all -- but it also breaks the three Cl3 formulas that divide by the
pseudoscalar or by a norm. Those are called out below.

Every grade is stored in **lexicographic** order, matching Cl3's `(xy, xz, yz)`:

| Type | Grade | Members | Geometry | Alias |
| --- | --- | --- | --- | --- |
| `Scalar` | 0 | `value` | | |
| `Vector` | 1 | `e0, e1, e2, e3` | a plane | `Plane` |
| `Bivector` | 2 | `e01, e02, e03, e12, e13, e23` | a line | `Line`, `Twist` |
| `Trivector` | 3 | `e012, e013, e023, e123` | a point | `Point` |
| `Quadvector` | 4 | `e0123` | the pseudoscalar | |
| `Multivector` | mixed | all sixteen components | | |
| `Rotor` | 0 + 2 (Euclidean) | `scalar, e12, e13, e23` | rotation about an origin axis | |
| `Translator` | 0 + 2 (ideal) | `scalar, e01, e02, e03` | translation | |
| `Motor` | 0 + 2 + 4 | `scalar, bivector, quadvector` | any rigid motion | |

`Rotor` and `Translator` store raw components rather than a `Bivector<T>`, so
"the ideal part is zero" (or "the Euclidean part is zero") is structural rather
than a precondition a constructor could quietly violate. `grade2()` hands the
half back as a `Bivector` when one is needed.

### Euclidean and ideal

The components split in two everywhere:

- **Euclidean**: the ones without `e0` -- `e1,e2,e3`, `e12,e13,e23`, `e123`.
  These carry direction and size.
- **Ideal**: the ones with `e0`. These carry position.

An object whose Euclidean part is zero is *ideal*: it lies at infinity. The
ideal plane is `e0`, an ideal point is a pure direction, an ideal line is a line
at infinity. `is_ideal()` tests for it, and `norm` and `ideal_norm` measure the
two halves separately.

### Products

Of the 256 entries in the 16×16 Cayley table, **64 vanish** because the product
contains `e0` twice; the remaining 192 take the sign of the permutation that
sorts the concatenated indices. `tests/test_pga.cpp` checks all 256 against an
independent implementation of that rule, so a hand-edited term in the table
fails the build.

Worth knowing:

```
e01 e12 = +e02        e12 e01 = -e02        e123 e123 = -1
e123 e0123 = +e0      e0123 e123 = -e0      e0123^2 = 0
e012 e123 = -e03      e012 e013 = 0         e0 e0 = 0
```

`e0123` commutes with every bivector and anticommutes with every vector, and

```
e12 e0123 = -e03      e13 e0123 = +e02      e23 e0123 = -e01
```

which is the coupling that `exp`, `log` and `screw` are built on.

**Grade of a product:** `Vector · Vector` → 0 + 2, `Bivector · Bivector` →
0 + 2 + 4, and **`Trivector · Trivector` → 0 + 2**. That last one differs from
Cl3, where two pseudoscalars multiply to a pure scalar. Here two points share
only two basis indices, so a bivector survives; `Trivector * Trivector` returns
a `Multivector`.

## Geometric primitives

Use the factories in `primitives.hpp` rather than writing components by hand --
the signs are not memorable, and that is the point of the factories existing.

```cpp
ga::plane(a, b, c, d)        // a x + b y + c z + d = 0  ->  d e0 + a e1 + b e2 + c e3
ga::point(x, y, z)           // e123 - x e023 + y e013 - z e012
ga::ideal_point(dx, dy, dz)  // the same, with zero weight: a direction
ga::origin<double>()         // e123
ga::ideal_plane<double>()    // e0
ga::line_through_origin(dx, dy, dz)   // dz e12 - dy e13 + dx e23
```

The x, y and z axes are therefore `e23`, `-e13` and `e12`.

A point's `e123` component is its **weight**; `position()` divides it out.
`point(x,y,z)` is the triple meet of the planes `X = x`, `Y = y`, `Z = z`, which
is what makes the incidence test exact:

```
plane ^ point = (a x + b y + c z + d) e0123
```

zero exactly when the point lies on the plane. (`point ^ plane` is its negation.)

A line's Euclidean part `(e12, e13, e23)` is its **direction** and its ideal
part `(e01, e02, e03)` is its **moment** -- `direction()` and `moment()` read
them back. The z-parallel line through `(1,0,0)` is `e12 - e02`.

### Meet and join

The wedge is the **meet**: two planes wedge to the line they share, a plane and
a line to the point they share, and `plane ^ point` to the incidence scalar
above. The dual operation -- the object *spanned* by two others -- is the
**regressive product**, `join()`:

```
join(a, b) = undual(dual(a) ^ dual(b))
```

`join(point, point)` is the line through both, whose Euclidean norm is their
distance; `join(point, line)` is the plane through both, and it vanishes when
the point lies on the line. Both are named functions, never operators: `^` is
the meet, and inventing a `&` for the join would read as "and" to most C++
programmers.

## Duality is a complement, not a product

**This is the first thing `e0^2 = 0` breaks.** In Cl3, `dual(A) = A e123`. Here
`e0123` squares to zero, so multiplying by it destroys information and cannot be
undone. PGA uses the **right complement** instead: `dual(a)` is the blade `b`
with `a ^ b = +e0123`.

| blade | dual | | blade | dual |
| --- | --- | --- | --- | --- |
| `1` | `+e0123` | | `e12` | `+e03` |
| `e0` | `+e123` | | `e13` | `-e02` |
| `e1` | `-e023` | | `e23` | `+e01` |
| `e2` | `+e013` | | `e012` | `+e3` |
| `e3` | `-e012` | | `e013` | `-e2` |
| `e01` | `+e23` | | `e023` | `+e1` |
| `e02` | `-e13` | | `e123` | `-e0` |
| `e03` | `+e12` | | `e0123` | `+1` |

`dual(dual(x))` is `+x` on grades 0, 2 and 4 and `-x` on grades 1 and 3, so
`undual()` -- the exact inverse -- is `dual` on even grades and `-dual` on odd
ones. It swaps the Euclidean and ideal halves, which is why `ideal_norm(x)`
equals `norm(dual(x))`.

## Norms

`norm(x)` is `sqrt(<x ~x>_0)`: the root of the sum of squares of the
**non-`e0`** components. It is the Euclidean size -- the length of a plane's
normal, of a line's direction, a point's weight -- and it is **zero for every
ideal element**, including the pseudoscalar. `ideal_norm(x)` measures the other
half. `squared_norm` skips the square root.

`normalize()` divides by `norm` and asserts on a zero one (the Cl3 shape).
There is no `normalize(Quadvector)`: its norm is always zero.

## Inverses

**The second thing `e0^2 = 0` breaks.** Not everything is invertible: `e0`,
`e0123` and every ideal element square to zero. The library returns a zero
object rather than a NaN in those cases, as Cl3's multivector inverse does.

- `Vector`, `Trivector`, `Rotor`, `Translator`: `x ~x` is a scalar, so
  `inverse(x) = reverse(x) / squared_norm(x)`, exactly as in Cl3.
- `Bivector`, `Motor`: `x ~x` is a **study number** `a + b e0123`, not a scalar.
  Since `e0123^2 = 0` it inverts as `(a - b e0123) / a^2`, so
  `inverse(x) = reverse(x) (a - b e0123) / a^2`. **Dropping the `b` term gives a
  wrong answer for any motor that is not already normalised**, which is the
  easiest mistake to make here.
- `Multivector`: the Hitzer–Sangwine formula for n = 4. With
  `N = m conjugate(m)` (grades 0, 3 and 4 only) and `N'` the same with grades 3
  and 4 negated, `N N'` is a pure scalar and
  `inverse(m) = conjugate(m) N' / <N N'>_0`. The derivation uses only the
  grading, so it survives the degenerate metric.
- `Quadvector`: no inverse; the overload does not exist.

Normalising a motor has the same wrinkle: `normalize(Motor)` divides by the
square root of the whole study number, `sqrt(a) (1 + (b/2a) e0123)`, so that
`m ~m = 1` exactly. Dividing by `norm(m)` alone leaves an `e0123` residue.

## Involutions

Metric-independent, so these are Cl3's signs extended to grade 4:

| | 0 | 1 | 2 | 3 | 4 |
| --- | --- | --- | --- | --- | --- |
| `reverse` | + | + | − | − | + |
| `involute` | + | − | + | − | + |
| `conjugate` | + | − | − | + | + |

## Motions

`sandwich(x, m)` applies `M x ~M` -- **object first, versor second**, the same
argument order as Cl3. One formula moves planes, lines, points and whole
multivectors, which is the practical payoff of PGA. `rotate`, `translate` and
`transform` are the same function under friendlier names.

### Translators

```
translator(dx, dy, dz) = 1 - (dx e01 + dy e02 + dz e03) / 2
```

The half is what the two-sided sandwich needs, the same as a rotor's half
angle. Normalised translators compose by **adding** their ideal parts, because
every ideal bivector squares to zero.

### Rotors and motors

`exp(-(θ/2) L)` rotates by `θ` about the unit line `L`. That is the ordinary
right-hand rule: `+90°` about `+z` takes `point(1,0,0)` to `point(0,1,0)` and
the plane `x = 0` to the plane `y = 0` -- identical to Cl3's
`rotor_from_axis_angle`.

What PGA adds is that **`L` need not pass through the origin**:
`motor_from_line_angle(L, θ)` rotates about any line, and
`screw(L, θ, d)` also translates `d` along it. Chasles' theorem says every
rigid motion is one of those, so every motor is a screw.

`motor_between(p, q) = normalize(1 + q p)` takes one plane onto another -- the
same half-angle construction as Cl3's `rotor_between`, with a half turn about a
non-degenerate line chosen when the planes are opposite.

### exp and log

`exp` does **not** collapse to `cos + sin` as it does in Cl3, because a general
bivector does not square to a scalar: `B^2 = -θ^2 + 2 t e0123`. Since
`e0123^2 = 0` the study numbers behave like dual numbers, and each of `cos` and
`sin` picks up a derivative term. With `B_e` the Euclidean part and
`B_e I` mapping `(e12, e13, e23)` onto `(-e03, +e02, -e01)`:

```
θ  = |B_e|                              (zero => exp(B) = 1 + B exactly)
b4 = 2 (e01 e23 - e02 e13 + e03 e12)    (= <B^2>_4)
d  = -b4 / (2θ)                         (the pitch)
s  = sin(θ)/θ,  c = cos(θ)

exp(B) = c + s B + (d/θ)(c - s) (B_e I) - d sin(θ) e0123
```

and its inverse, for a unit motor with scalar `s0`, bivector `B`, `e0123` part `p`:

```
n = |B_e|                               (zero => log(m) = B / s0, a translator)
θ = atan2(n, s0),  d = -p / n

log(m) = (θ/n) B + (d/n)(1 - θ cos(θ)/n) (B_e I)
```

Both were checked against a forty-term power series. `exp(B)` is always a unit
motor. The round trip is exact inside the principal range and wraps outside it,
as Cl3's does.

`slerp(a, b, t) = a exp(t log(~a b))`, negating `b` first if the two point away
from each other. Interpolating this way moves along a **single screw**: the
rotation and translation stay coupled, so the path of a point is a helix rather
than a corner.

## Geometry

```cpp
reflect(x, plane)                  // planes and lines flip sign, points do not
project(point, plane)              // and (point, line), (line, plane)
distance(point, point)             // norm(join) / the weights
distance(point, plane)             // SIGNED: which side you are on
distance(point, line)              // norm(join) / weight / line norm
angle(plane, plane)                // acos of the normalized dot
angle(line, line)                  // acos of minus the normalized scalar product
```

Projection is `(x ⌊ y) inverse(y)`, the same shape as Cl3's. Reflection needs a
grade-dependent sign -- planes and lines take `-p x p⁻¹`, points take
`+p x p⁻¹` -- which is what keeps a reflected point's weight positive.

The offsets do all the work you would expect: two parallel planes at different
offsets are at zero degrees to each other but reflect a point to different
places, and a line that misses the origin is a perfectly ordinary argument
everywhere.

`distance` and `angle` return a raw `T` rather than a `Scalar<T>` and are
`inline` rather than `constexpr`, since they use `std::acos` and `std::sqrt` --
the same exception `approx_equal` makes.

## What differs from Cl3, in one list

1. **`dual` is the complement**, not a product with the pseudoscalar, because
   `e0123^2 = 0`. `dual(dual(x))` flips sign on grades 1 and 3.
2. **`Trivector * Trivector` is a `Multivector`**, grades 0 and 2 -- not a
   `Scalar`.
3. **`norm` ignores every `e0` component**, so ideal elements have zero norm and
   `ideal_norm` exists to measure them.
4. **A motor's inverse and normalisation need the study number**; the naive
   `reverse(m)/|m|²` is wrong for a motor that has drifted.
5. **`Quadvector` has no inverse** and no `normalize`.
6. **`exp` carries derivative terms**; only the pure-Euclidean case reduces to
   Cl3's `cos + sin B/|B|`.
7. **Grade 1 is a plane, grade 3 a point.** `Plane`, `Line`, `Point` and `Twist`
   are aliases that say so.

## Known gaps

- No `exp`/`log` for a general multivector, only for the bivector/motor pair.
  Same deliberate choice as Cl3.
- No contractions involving `Multivector`, `Rotor`, `Translator` or `Motor`;
  the family covers the blade types, as in Cl3.
- `slerp` inherits `log`'s wrapping outside the principal range.
- `normalize(Bivector)` divides by the Euclidean norm only -- a "unit line"
  means unit *direction*, and the moment scales along with it.
