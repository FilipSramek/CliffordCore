#pragma once

#include <cmath>
#include <limits>
#include "../scalar.hpp"
#include "../vector3.hpp"
#include "../bivector3.hpp"
#include "../rotor3.hpp"
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

namespace CliffordCore
{
    template<typename T>
    /**
     * @brief The rotor that represents no rotation.
     * @return The identity rotor, scalar 1 with a zero bivector.
     */
    constexpr Rotor3<T> identity_rotor() {
        return Rotor3<T>(Scalar<T>(1), Bivector3<T>(0, 0, 0));
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
    constexpr Rotor3<T> rotor_from_axis_angle(const Vector3<T>& axis, T angle) {
        if (norm(axis).value == T(0)) {
            return identity_rotor<T>();
        }
        // dual() is a signed permutation, so the dual of a unit vector is
        // already a unit bivector -- the plane perpendicular to the axis.
        const Bivector3<T> plane = dual(normalize(axis));
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
    constexpr Rotor3<T> rotor_between(const Vector3<T>& from, const Vector3<T>& to) {
        if (norm(from).value == T(0) || norm(to).value == T(0)) {
            return identity_rotor<T>();
        }

        const Vector3<T> a = normalize(from);
        const Vector3<T> b = normalize(to);
        const T alignment = (a | b).value;

        // Opposite directions: every plane containing `a` is a valid half turn,
        // so pick one that is definitely not degenerate.
        if (alignment < T(-1) + std::numeric_limits<T>::epsilon() * T(8)) {
            const Vector3<T> reference =
                (std::abs(a.x) < T(0.9)) ? Vector3<T>(1, 0, 0) : Vector3<T>(0, 1, 0);
            // A unit bivector with zero scalar part is exactly a half turn.
            return Rotor3<T>(Scalar<T>(0), normalize(a ^ reference));
        }

        const Multivector3<T> product = geometric_product(b, a);
        return normalize(Rotor3<T>(product.scalar + Scalar<T>(1), product.bivector));
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
    constexpr Rotor3<T> slerp(const Rotor3<T>& from, const Rotor3<T>& to, T t) {
        const Rotor3<T> a = normalize(from);
        Rotor3<T> b = normalize(to);

        const T alignment = a.scalar.value * b.scalar.value
                          + a.bivector.xy * b.bivector.xy
                          + a.bivector.xz * b.bivector.xz
                          + a.bivector.yz * b.bivector.yz;
        if (alignment < T(0)) {
            b = Rotor3<T>(Scalar<T>(-b.scalar.value),
                          Bivector3<T>(-b.bivector.xy, -b.bivector.xz, -b.bivector.yz));
        }

        // reverse(a) is the inverse of a unit rotor, so this is the relative
        // rotation from a to b, scaled by t and reapplied.
        const Bivector3<T> between = log(reverse(a) * b);
        return a * exp(between * Scalar<T>(t));
    }
} // namespace CliffordCore
