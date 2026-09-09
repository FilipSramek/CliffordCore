# Rotations

Rotors are the reason to reach for geometric algebra in 3D. They compose by
multiplication, never gimbal-lock, interpolate cleanly, and carry no redundant
degrees of freedom beyond a single normalisation constraint.

A rotor is a scalar plus a bivector — `Rotor3<T>` holds exactly those two parts.
Concretely it is the product of two unit vectors, which is the fact that explains
everything else on this page.

Runnable companions: `examples/03_rotations.cpp` and
`examples/04_interpolation.cpp`.

## Building a rotor

| From | Call |
| --- | --- |
| nothing (no rotation) | `identity_rotor<double>()` |
| an axis and an angle | `rotor_from_axis_angle(axis, radians)` |
| "point this at that" | `rotor_between(from, to)` |
| a bivector (a plane and an angle) | `exp(bivector)` |
| two unit vectors | `rotor_product(a, b)` |

`rotor_from_axis_angle` normalises the axis internally, so any non-zero length
works, and returns the identity for a zero axis. `rotor_between` returns the
identity for parallel inputs and a genuine half turn for antiparallel ones —
that case is ambiguous (every plane containing the vector is a valid half turn)
so the library picks one rather than returning a degenerate rotor.

## The sign that bites you

These two are both quarter turns in the xy plane, **and they go opposite ways**:

```cpp
auto viaExp  = exp(Bivector3<double>(kPi / 4, 0, 0));
auto viaAxis = rotor_from_axis_angle(Vector3<double>(0, 0, 1), kPi / 2);

rotate(e1, viaExp);    // (0, -1, 0)
rotate(e1, viaAxis);   // (0,  1, 0)
```

Three separate things are in play:

1. **`exp(theta * B)` rotates by `2 * theta`.** A quarter turn is `exp(pi/4 * B)`.
   The half-angle is the same one that appears in quaternions, and for the same
   reason: the rotor is applied twice, once on each side.
2. **`e12` spans the *xy plane*, so it rotates about z — not about x.** Pairing
   the `xy` component with the x axis is the classic mistake. The axis of a
   bivector is its dual, and `dual(e12) = -e3`.
3. **`rotor_from_axis_angle` is `exp(-(angle/2) * dual(axis))`.** That minus sign
   is precisely what makes it come out as the ordinary right-hand rule.

If you build rotors with `rotor_from_axis_angle` and apply them with `rotate`,
everything behaves the way you expect and you never have to think about this.
It only surfaces when you hand-build a bivector and call `exp` on it.

## Applying a rotor

```cpp
Vector3<double> turned = rotate(v, r);      // or sandwich(v, r)
```

**The vector comes first, the rotor second.** `rotate` is a friendlier alias for
`sandwich`, which computes `R v reverse(R)`.

Spelled out longhand, that sandwich is:

```cpp
Vector3<double> turned = grade1(r * v * reverse(r));
```

which exercises `Rotor3 * Vector3` and `Multivector3 * Rotor3` from
`mixed_products.hpp`. The scalar, bivector and trivector parts of that product
are all zero — a vector sandwiched by a rotor is always another vector.
`sandwich()` is the closed form of exactly this expression, which is why it
exists: same answer, far less arithmetic.

## Composing

```cpp
Rotor3<double> both = second * first;   // apply `first`, then `second`
```

`Rotor3 * Rotor3` stays a `Rotor3` — composition never widens. And rotations do
not commute, so neither does the product: `Rz * Rx` and `Rx * Rz` send `e1` to
different places.

## Undoing

For a **unit** rotor, the inverse is just the reverse:

```cpp
reverse(r)     // flip the bivector part
inverse(r)     // reverse(r) / squared_norm(r)
```

`inverse` additionally divides by the squared norm, which only matters once the
rotor has drifted off unit length.

## Drift

Every composition accumulates a little floating-point error and the rotor slowly
stops being unit length — at which point it scales as well as rotates. Measured
on `|norm - 1|` after repeated products (see `examples/03_rotations.cpp`):

| `T` | after 1 000 | after 1 000 000 |
| --- | --- | --- |
| `double` | 8.9e-16 | 2.6e-12 |
| `float` | 1.8e-05 | 1.7e-02 |

In `double` this is negligible for almost any workload. In `float` it is a
visible 1.7% scaling error after a million compositions — a unit vector comes
back with length 0.9665. The fix is `normalize(r)`, applied periodically.

## Interpolating

```cpp
Rotor3<double> midway = slerp(from, to, 0.5);
```

Endpoints are exact and the intermediate steps are evenly spaced in angle. A
rotor and its negation represent the *same* rotation, so `slerp` negates its
second operand when the two point away from each other — without that you would
take the 270-degree path instead of the 90-degree one.

## Recovering an axis and angle

```cpp
Bivector3<double> b = log(r);
double angle        = 2.0 * norm(b).value;
Vector3<double> axis = dual(normalize(b));
```

The `2` is the half-angle from `exp` reappearing. The `dual` undoes the `dual`
that `rotor_from_axis_angle` applied — up to sign, and the signs work out
because `dual(dual(x)) == -x` cancels against the `-(angle/2)`.

`log` clamps the rotor's scalar part into `[-1, 1]` before calling `std::acos`.
That matters: composition drifts, `acos` is undefined outside that range, and
without the clamp a drifted rotor would yield `NaN`.

## Why the sandwich has that shape

Reflecting a vector in the plane with unit normal `n` is `-n v n`. Reflect in
`n`, then in `m`:

```text
m (-n v n) m  =  (m n) v (n m)  =  R v reverse(R)      where R = m n
```

So a rotor **is** the product of two unit vectors, and the sandwich **is** two
reflections in a row. The rotation angle is twice the angle between the two
planes. `examples/05_geometry_toolkit.cpp` prints both sides of this and shows
they match.

That is also why `rotor_product(a, b)` exists as a named function: it is the
geometric product of two vectors, packed into the type that represents the
rotation they generate.
