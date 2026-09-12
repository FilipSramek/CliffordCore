#pragma once

/**
 * @file cliffordcore/cl3/operations/sandwich.hpp
 * @brief Applying a rotor to a vector via the sandwich product.
 *
 * sandwich(x, r) applies R x ~R, with the object first and the rotor
 * second; rotate is the same function under a friendlier name. Overloads
 * cover vectors, bivectors, trivectors and multivectors -- the trivector
 * case is the identity, since the pseudoscalar is central.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <type_traits>
#include "../vector.hpp"
#include "../rotor.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../multivector.hpp"
#include "grade.hpp"
#include "mixed_products.hpp"
#include "reverse.hpp"

namespace CliffordCore::Cl3
{
    template<typename T>
    /**
     * @brief Computes the sandwich product of a 3D vector and a 3D rotor (R*v*Reverse(R)).
     * @param v The vector to be transformed.
     * @param r The rotor to perform the transformation with.
     * @return The resulting vector after the sandwich product.
     */
    constexpr Vector<T> sandwich(const Vector<T>& v, const Rotor<T>& r) {
        const T s  = r.scalar.value;
        const T xy = r.bivector.xy;
        const T xz = r.bivector.xz;
        const T yz = r.bivector.yz;

        // Expanding R v ~R. Note which axis each bivector turns: e12 (xy) spans
        // the xy plane, so it rotates about z -- not about x. Pairing xy with
        // the x axis is the easy mistake here.
        return Vector<T>(
            (s*s - xy*xy - xz*xz + yz*yz) * v.x
                + 2 * (s*xy - xz*yz) * v.y
                + 2 * (s*xz + xy*yz) * v.z,

            2 * (-s*xy - xz*yz) * v.x
                + (s*s - xy*xy + xz*xz - yz*yz) * v.y
                + 2 * (s*yz - xy*xz) * v.z,

            2 * (-s*xz + xy*yz) * v.x
                + 2 * (-s*yz - xy*xz) * v.y
                + (s*s + xy*xy - xz*xz - yz*yz) * v.z
        );
    }

    /**
     * @brief Rotates a 3D vector using a 3D rotor.
     * @param v The vector to be rotated.
     * @param r The rotor to perform the rotation with.
     * @return The resulting vector after rotation.
     */
    template<typename T>
    constexpr Vector<T> rotate(const Vector<T>& v, const Rotor<T>& r) {
        return sandwich(v, r);
    }

    template<typename T>
    /**
     * @brief Rotates a bivector by a rotor.
     * @param b The bivector to rotate.
     * @param r The rotor to rotate with.
     * @return The rotated bivector, R b reverse(R).
     *
     * Rotating a plane rather than a direction. Needed whenever the thing being
     * turned is an orientation, an angular velocity, or a surface element.
     */
    constexpr Bivector<T> sandwich(const Bivector<T>& b, const Rotor<T>& r) {
        return grade2(r * b * reverse(r));
    }

    template<typename T>
    /**
     * @brief Rotates a trivector by a rotor.
     * @param t The trivector to rotate.
     * @param r The rotor to rotate with.
     * @return The trivector, scaled by the rotor's squared norm. The
     *         pseudoscalar commutes with everything in 3D, so a unit rotor
     *         leaves it exactly unchanged -- volume has no orientation to turn.
     */
    constexpr Trivector<T> sandwich(const Trivector<T>& t, const Rotor<T>& r) {
        return grade3(r * t * reverse(r));
    }

    template<typename T>
    /**
     * @brief Rotates a multivector by a rotor.
     * @param m The multivector to rotate.
     * @param r The rotor to rotate with.
     * @return The rotated multivector. Every grade is carried along, so no
     *         projection is needed here.
     */
    constexpr Multivector<T> sandwich(const Multivector<T>& m, const Rotor<T>& r) {
        return r * m * reverse(r);
    }

    template<typename T>
    /**
     * @brief Rotates a bivector by a rotor.
     * @param b The bivector to rotate.
     * @param r The rotor to rotate with.
     * @return The rotated bivector.
     */
    constexpr Bivector<T> rotate(const Bivector<T>& b, const Rotor<T>& r) {
        return sandwich(b, r);
    }

    template<typename T>
    /**
     * @brief Rotates a trivector by a rotor.
     * @param t The trivector to rotate.
     * @param r The rotor to rotate with.
     * @return The rotated trivector.
     */
    constexpr Trivector<T> rotate(const Trivector<T>& t, const Rotor<T>& r) {
        return sandwich(t, r);
    }

    template<typename T>
    /**
     * @brief Rotates a multivector by a rotor.
     * @param m The multivector to rotate.
     * @param r The rotor to rotate with.
     * @return The rotated multivector.
     */
    constexpr Multivector<T> rotate(const Multivector<T>& m, const Rotor<T>& r) {
        return sandwich(m, r);
    }
} // namespace CliffordCore::Cl3
