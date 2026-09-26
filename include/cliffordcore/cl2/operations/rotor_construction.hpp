#pragma once

/**
 * @file cliffordcore/cl2/operations/rotor_construction.hpp
 * @brief Building rotors from angles and vector pairs, reading their angle
 *        back, and interpolating between them.
 *
 * rotor_from_angle is counter-clockwise for a positive angle -- +90 degrees
 * takes e1 to e2, exactly as Cl(3,0) does about +z. rotor_between builds the
 * rotor carrying one vector to another, rotor_angle reads the angle back
 * out, and slerp interpolates along the shortest arc.
 *
 * @author Filip Sramek
 * @version 0.3.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <cmath>
#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../rotor.hpp"
#include "dot_product.hpp"
#include "exp.hpp"
#include "log.hpp"
#include "geometric_product.hpp"
#include "norm.hpp"
#include "normalize.hpp"
#include "reverse.hpp"
#include "wedge_product.hpp"

// Ways to build a rotor without assembling a bivector by hand.
//
// There is no axis in two dimensions -- every rotation happens in the one
// plane e12 -- so where Cl(3,0) has rotor_from_axis_angle(axis, angle), this
// algebra has rotor_from_angle(angle). It is the one construction function
// whose name does not survive a namespace switch, deliberately: an axis
// argument here would be a fiction.
//
// Orientation convention: exp(theta * e12) rotates by -2*theta, and sandwich()
// applies it as R v reverse(R). A counter-clockwise rotation by `angle` is
// therefore exp(-(angle/2) * e12), the same sign Cl(3,0) uses about +z.

namespace CliffordCore::Cl2
{
    template<typename T>
    /**
     * @brief The rotor that represents no rotation.
     * @return The identity rotor, scalar 1 with a zero bivector.
     */
    constexpr Rotor<T> identity_rotor() {
        return Rotor<T>(Scalar<T>(1), Bivector<T>(0));
    }

    template<typename T>
    /**
     * @brief Builds a rotor from an angle.
     * @param angle The rotation angle in radians, counter-clockwise for a
     *        positive value: +pi/2 takes e1 to e2.
     * @return The rotor performing that rotation, cos(angle/2) - sin(angle/2) e12.
     */
    constexpr Rotor<T> rotor_from_angle(T angle) {
        return exp(Bivector<T>(-angle / T(2)));
    }

    template<typename T>
    /**
     * @brief Reads back the angle a rotor rotates by.
     * @param r The rotor. Need not be unit length.
     * @return The rotation angle in radians, in (-pi, pi], counter-clockwise
     *         positive. A rotor and its negation are the same rotation and give
     *         the same angle, so rotor_angle(rotor_from_angle(a)) == a only for
     *         a in that range; outside it, the equivalent angle comes back.
     *
     * Computed as the angle e1 is turned to, which is atan2(-2sb, s^2 - b^2) for
     * r = s + b e12. That form is invariant under r -> -r and under scaling, so
     * it needs neither a normalize() nor a wrap into range.
     */
    inline T rotor_angle(const Rotor<T>& r) {
        const T s = r.scalar.value;
        const T b = r.bivector.xy;
        // + T(0) turns a -0.0 into +0.0, so an exact half turn reports +pi
        // rather than -pi. See the note in geometry.hpp.
        return std::atan2(T(-2) * s * b + T(0), s * s - b * b);
    }

    template<typename T>
    /**
     * @brief Builds the rotor that rotates one vector onto another.
     * @param from The starting direction. Normalized internally.
     * @param to The target direction. Normalized internally.
     * @return A unit rotor r with sandwich(from, r) parallel to to.
     *
     * Built as normalize(1 + to*from), the half-angle construction, evaluated
     * so that it stays accurate to a few epsilons however close the two
     * directions are to opposite. When they are exactly opposite that
     * expression is zero. There is no choice
     * of plane to make in two dimensions -- only which way round -- so the
     * counter-clockwise half turn, rotor_from_angle(pi), is returned.
     */
    constexpr Rotor<T> rotor_between(const Vector<T>& from, const Vector<T>& to) {
        if (norm(from).value == T(0) || norm(to).value == T(0)) {
            return identity_rotor<T>();
        }

        const Vector<T> a = normalize(from);
        const Vector<T> b = normalize(to);
        const T d = (a | b).value;

        // The rotor is normalize((1 + a.b) + b ^ a), but both parts are
        // computed so that neither cancels when a and b are nearly opposite --
        // the naive 1 + a.b loses a digit for every factor of ten the angle
        // gets closer to pi, which float runs out of at about 1e-3 radians.
        //
        //   b ^ a == (a + b) ^ a exactly, since a ^ a = 0, and a + b is
        //   computed without rounding when b is close to -a.
        //   1 + a.b == |a ^ b|^2 / (1 - a.b) exactly, since
        //   |a ^ b|^2 = 1 - (a.b)^2 for unit vectors; used only when a.b < 0,
        //   where 1 - a.b is close to 2 and nothing cancels.
        const Bivector<T> plane = (a + b) ^ a;
        const T s = d >= T(0) ? T(1) + d : (plane.xy * plane.xy) / (T(1) - d);

        // Exactly opposite: both parts vanish, and a half turn either way is
        // correct. Take the counter-clockwise one, written exactly.
        if (s == T(0) && plane.xy == T(0)) {
            return Rotor<T>(Scalar<T>(0), Bivector<T>(-1));
        }

        return normalize(Rotor<T>(Scalar<T>(s), plane));
    }

    template<typename T>
    /**
     * @brief Interpolates between two rotors along the shortest arc.
     * @param from The rotor at t = 0.
     * @param to The rotor at t = 1.
     * @param t The interpolation parameter, normally in [0, 1].
     * @return The interpolated rotor, r0 * exp(t * log(reverse(r0) * r1)).
     *
     * A rotor and its negation describe the same rotation, so if the two rotors
     * point away from each other the second is negated first and the result
     * takes the short way round.
     */
    constexpr Rotor<T> slerp(const Rotor<T>& from, const Rotor<T>& to, T t) {
        const Rotor<T> a = normalize(from);
        Rotor<T> b = normalize(to);

        const T alignment = a.scalar.value * b.scalar.value
                          + a.bivector.xy * b.bivector.xy;
        if (alignment < T(0)) {
            b = Rotor<T>(Scalar<T>(-b.scalar.value), Bivector<T>(-b.bivector.xy));
        }

        // reverse(a) is the inverse of a unit rotor, so this is the relative
        // rotation from a to b, scaled by t and reapplied.
        const Bivector<T> between = log(reverse(a) * b);
        return a * exp(between * Scalar<T>(t));
    }
} // namespace CliffordCore::Cl2
