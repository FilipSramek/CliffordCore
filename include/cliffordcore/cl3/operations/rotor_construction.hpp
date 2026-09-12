#pragma once

/**
 * @file cliffordcore/cl3/operations/rotor_construction.hpp
 * @brief Building rotors from axes, angles, vector pairs, and interpolation.
 *
 * rotor_from_axis_angle follows the ordinary right-hand rule -- +90 degrees
 * about +z takes e1 to e2. rotor_between builds the rotor carrying one
 * vector to another, and slerp interpolates along the shortest arc.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <cmath>
#include <limits>
#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../rotor.hpp"
#include "dual.hpp"
#include "exp.hpp"
#include "log.hpp"
#include "geometric_product.hpp"
#include "norm.hpp"
#include "normalize.hpp"
#include "reverse.hpp"
#include "wedge_product.hpp"

// Ways to build a rotor without assembling a bivector by hand.
//
// Orientation convention: exp(theta * B) rotates by 2*theta in the plane B, and
// sandwich() applies it as R v reverse(R). Working that through, a rotation of
// `angle` about `axis` is exp(-(angle/2) * dual(axis)), which turns out to be
// the ordinary right-hand rule: +90 degrees about +z takes e1 to e2.

namespace CliffordCore::Cl3
{
    template<typename T>
    /**
     * @brief The rotor that represents no rotation.
     * @return The identity rotor, scalar 1 with a zero bivector.
     */
    constexpr Rotor<T> identity_rotor() {
        return Rotor<T>(Scalar<T>(1), Bivector<T>(0, 0, 0));
    }

    template<typename T>
    /**
     * @brief Builds a rotor from an axis and an angle.
     * @param axis The axis to rotate about. It is normalized internally, so any
     *        non-zero length works.
     * @param angle The rotation angle in radians, right-handed about the axis.
     * @return The rotor performing that rotation, or the identity if the axis
     *         has zero length.
     */
    constexpr Rotor<T> rotor_from_axis_angle(const Vector<T>& axis, T angle) {
        if (norm(axis).value == T(0)) {
            return identity_rotor<T>();
        }
        // dual() is a signed permutation, so the dual of a unit vector is
        // already a unit bivector -- the plane perpendicular to the axis.
        const Bivector<T> plane = dual(normalize(axis));
        return exp(plane * Scalar<T>(-angle / T(2)));
    }

    template<typename T>
    /**
     * @brief Builds the rotor that rotates one vector onto another.
     * @param from The starting direction. Normalized internally.
     * @param to The target direction. Normalized internally.
     * @return A unit rotor r with sandwich(from, r) parallel to to.
     *
     * Built as normalize(1 + to*from), the half-angle construction. When the two
     * directions are opposite the rotation plane is ambiguous, so a perpendicular
     * plane containing `from` is chosen and the result is a half turn.
     */
    constexpr Rotor<T> rotor_between(const Vector<T>& from, const Vector<T>& to) {
        if (norm(from).value == T(0) || norm(to).value == T(0)) {
            return identity_rotor<T>();
        }

        const Vector<T> a = normalize(from);
        const Vector<T> b = normalize(to);
        const T alignment = (a | b).value;

        // Opposite directions: every plane containing `a` is a valid half turn,
        // so pick one that is definitely not degenerate.
        if (alignment < T(-1) + std::numeric_limits<T>::epsilon() * T(8)) {
            const Vector<T> reference =
                (std::abs(a.x) < T(0.9)) ? Vector<T>(1, 0, 0) : Vector<T>(0, 1, 0);
            // A unit bivector with zero scalar part is exactly a half turn.
            return Rotor<T>(Scalar<T>(0), normalize(a ^ reference));
        }

        const Multivector<T> product = geometric_product(b, a);
        return normalize(Rotor<T>(product.scalar + Scalar<T>(1), product.bivector));
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
                          + a.bivector.xy * b.bivector.xy
                          + a.bivector.xz * b.bivector.xz
                          + a.bivector.yz * b.bivector.yz;
        if (alignment < T(0)) {
            b = Rotor<T>(Scalar<T>(-b.scalar.value),
                          Bivector<T>(-b.bivector.xy, -b.bivector.xz, -b.bivector.yz));
        }

        // reverse(a) is the inverse of a unit rotor, so this is the relative
        // rotation from a to b, scaled by t and reapplied.
        const Bivector<T> between = log(reverse(a) * b);
        return a * exp(between * Scalar<T>(t));
    }
} // namespace CliffordCore::Cl3
