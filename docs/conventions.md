# Conventions

Every sign in CliffordCore follows from a handful of choices made once. This page
collects them. If a result's sign looks wrong, the answer is almost always here.

## The basis

Vectors are `e1 e2 e3` with `e_i^2 = +1`. Bivectors are stored as

```text
(xy, xz, yz)  =  (e1e2, e1e3, e2e3)
```

**Note `xz`, not the more common `zx`.** That single choice determines every
sign in the wedge product, the dual, and the rotation formulas. It is the first
thing to check when something comes out negated.

The pseudoscalar is `e123 = e1e2e3`. From `e_i^2 = +1` it follows that

```text
e12^2 = e13^2 = e23^2 = e123^2 = -1
```

Those four minus signs drive every inverse and dual formula below.

## Component names

| Type | Grade | Members |
| --- | --- | --- |
| `Scalar<T>` | 0 | `value` |
| `Vector<T>` | 1 | `x, y, z` |
| `Bivector<T>` | 2 | `xy, xz, yz` |
| `Trivector<T>` | 3 | `e123` |
| `Multivector<T>` | mixed | `scalar, vector, bivector, trivector` |
| `Rotor<T>` | 0 + 2 | `scalar, bivector` |

`Trivector`'s member is **`e123`**, not `value`. Reaching for `.value` on a
trivector is a recurring mistake.

## Return-type rules

**`operator*` always returns the most general type the product can produce.**
There is therefore exactly one `operator*` per operand pair. This is forced:
C++ cannot overload on return type, so a second `operator*(Vector, Vector)`
returning `Rotor` would make every `v * v` ambiguous.

**Narrower spellings are named functions:**

| Instead of | Use | To get |
| --- | --- | --- |
| `a * b` | `rotor_product(a, b)` | `Rotor` from two vectors |
| `s + b` | `rotor_sum(s, b)` | `Rotor` from a scalar and a bivector |
| — | `to_rotor(m)` | `Rotor` from a multivector (drops grades 1 and 3) |
| — | `to_multivector(r)` | `Multivector` from a rotor |

**Mixed-grade `+` and `-` always widen to `Multivector`** — including
`Scalar + Bivector`, even though that pair is exactly a rotor. Uniformity keeps
chaining predictable. Only binary pairs are defined; `a + b + c` parses as
`((a + b) + c)`, so chains of any length follow from the pairs alone.
Same-grade sums keep their own type: `Vector + Vector` is a `Vector`.

Every mixed-grade product routes through the single 8x8 multiplication table in
`geometric_product.hpp`. There is exactly one such table in the library; when
adding a product, delegate to it rather than deriving a closed form.

## The Trivector exception

Two pseudoscalars multiply to a **scalar**, not a trivector:

```cpp
Trivector<double> t(6), u(3);
t * u        // Scalar<double>(-18)   because e123 * e123 = -1
t / u        // Scalar<double>(2)     the two minus signs cancel
t * 2.0      // Trivector<double>(12)   scaling stays a trivector
t * Scalar<double>(2)   // Trivector<double>(12)
```

So `t * u` and `t * 2.0` deliberately differ in return type. `Trivector` is the
only grade type with a single-argument converting constructor, which is why it
carries extra raw-`T` overloads: without them a bare `2.0` could reach either
`operator*(Trivector)` or `operator*(Scalar<T>)` by one user-defined conversion
each, and the call would be ambiguous.

## The inner product

**`operator|` is the left contraction.** For two vectors that is exactly the dot
product it has always been, so nothing about `a | b` has changed; the contraction
is simply what it generalises to across grades:

```text
A _| B  =  <A B>_(s-r)        for blades of grade r and s
```

Where the grade would go negative the result is zero, so `Bivector | Vector` is
not provided -- you almost always meant `Vector | Bivector`.

The left contraction is chosen over the symmetric "fat dot" deliberately. It is
the product the projection formulas are written with -- `project(v, plane)` is
`(v _| B) * inverse(B)` -- and it is the one that keeps working in a degenerate
metric, which will matter when Cl(3,0,1) lands. Committing to that reading now
means `|` will not quietly change meaning later.

`scalar_product(a, b)` is the grade 0 part of any product, and
`right_contraction` is the reverse-mirror of the left one.

## Duality

`dual(A) = A * e123`. Working the basis through:

```text
dual(1)   = e123        dual(e1)  =  e23
dual(e12) = -e3         dual(e2)  = -e13
dual(e13) =  e2         dual(e3)  =  e12
dual(e23) = -e1         dual(e123) = -1
```

Because `e123^2 = -1`, **applying the dual twice negates**:

```cpp
dual(dual(a)) == -a     // correct, not a defect
```

The dual relates the wedge to the cross product you already know:

```cpp
dual(b ^ a) == a x b    // the classical cross product
dual(a ^ b) == -(a x b)
```

A plane and its normal axis are dual descriptions of the same thing.

*(The other common convention, `A * inverse(e123)`, negates every result and
makes `dual(a ^ b)` the cross product directly. This library uses `A * e123`.)*

## The three involutions

They come as a set, and each scales grade `k` by a sign depending only on `k`:

| | grade 0 | grade 1 | grade 2 | grade 3 |
| --- | --- | --- | --- | --- |
| `reverse` | + | + | - | - |
| `involute` | + | - | + | - |
| `conjugate` | + | - | - | + |

So `conjugate == reverse` composed with `involute`, in either order -- an
identity the test suite asserts rather than trusting the table. All three are
defined for all six types.

## Inverses

**Inverses carry the minus signs the squares imply:**

| Input | Inverse | Why |
| --- | --- | --- |
| `Vector` | `v / \|v\|^2` | `v^2 = +\|v\|^2` |
| `Bivector` | `-b / \|b\|^2` | `b^2 = -\|b\|^2` |
| `Trivector` | `-t / \|t\|^2` | `e123^2 = -1` |
| `Rotor` | `reverse(r) / \|r\|^2` | `r * reverse(r) = \|r\|^2` |
| `Multivector` | via Clifford conjugation | see below |

For a **unit** rotor the inverse is simply the reverse — which is why rotations
undo by reversing.

`inverse(Multivector)` is **not** `m / |m|^2`; that only works for a single
grade (a unit bivector would give `b * b = -1`). Instead, `m * conjugate(m)`
collapses to a scalar plus a pseudoscalar, `a + b*e123`. Since `e123^2 = -1`
that behaves like a complex number and inverts as `(a - b*e123) / (a^2 + b^2)`.
A non-invertible multivector returns a zero multivector rather than dividing by
zero.

## Orientation

This is the section that saves the most time.

- **`sandwich(v, r)` applies `R v reverse(R)`.** The **vector comes first**, the
  rotor second. `rotate(v, r)` is the same function under a friendlier name.
- **`exp(theta * B)` rotates by `2 * theta`** in the plane `B`. A quarter turn is
  `exp(pi/4 * B)`, not `exp(pi/2 * B)`.
- **`e12` spans the *xy plane*, so it rotates about z — not about x.** Pairing
  `xy` with the x axis is the classic mistake, and was once a real bug here.
- **`rotor_from_axis_angle` is the ordinary right-hand rule**: +90 degrees about
  +z takes `e1` to `e2`. It is `exp(-(angle/2) * dual(axis))`, and that minus
  sign is exactly what makes the handedness come out conventional.

The consequence worth internalising: `exp(pi/4 * e12)` and
`rotor_from_axis_angle({0,0,1}, pi/2)` are both quarter turns in the xy plane
**and they go opposite ways**. See `examples/03_rotations.cpp`, which prints them
side by side.

## Known gaps

- **`operator+` and `operator-` are undefined for `Rotor`** with anything but
  another `Rotor`, plus `Multivector + Rotor` which widens losslessly. Use
  `to_multivector(r)` first. `operator*`, by contrast, is defined for every pair.
- **No contractions between a `Multivector` and anything.** The contraction
  family covers the blade types, which is where it is well defined.
- **No general `exp`/`log`** for arbitrary multivectors -- only
  `exp(Bivector)` and `log(Rotor)`. The series does not close cleanly once odd
  grades are involved, and `log` becomes multivalued.
- **No `meet` / `join`.** In Cl(3,0) everything passes through the origin, which
  makes them much less interesting than they are in a projective algebra.

### A conversion that used to bite

`Rotor`'s constructor from `Multivector` is **explicit**. While it was
implicit, `r + m` compiled, narrowed `m`, and discarded grades 1 and 3 silently,
while `m + r` kept everything. Narrowing now has to be spelled out, with the
constructor or with `to_rotor(m)`.
