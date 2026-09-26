#pragma once

/**
 * @file cliffordcore/cl2/operations/sandwich.hpp
 * @brief Applying a rotor via the sandwich product.
 *
 * sandwich(x, r) applies R x ~R, with the object first and the rotor
 * second; rotate is the same function under a friendlier name. Overloads
 * cover vectors, bivectors and multivectors -- the bivector case is the
 * identity for a unit rotor, since the pseudoscalar commutes with every
 * rotor.
 *
 * @author Filip Sramek
 * @version 0.3.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <type_traits>
#include "../vector.hpp"
#include "../rotor.hpp"
#include "../bivector.hpp"
#include "../multivector.hpp"
#include "geometric_product.hpp"
#include "grade.hpp"
#include "mixed_products.hpp"
#include "reverse.hpp"

// The two-sided sandwich collapses to a one-sided product in two dimensions.
//
// The pseudoscalar anticommutes with vectors, so for R = s + b e12,
//
//     v ~R = v (s - b e12) = s v + b e12 v = (s + b e12) v = R v
//
// and therefore R v ~R = R R v = R^2 v. A rotation of a vector in the plane is
// just multiplication on the left by R^2 = (s^2 - b^2) + 2sb e12 -- which, read
// as complex numbers, is the familiar "multiply by e^(i theta)". The identity
// holds for any rotor, unit or not; the tests assert it directly.
//
// sandwich() is still written and named as a sandwich, so code carries over
// from Cl(3,0) unchanged and reads the same in both algebras.

namespace CliffordCore::Cl2
{
    template<typename T>
    /**
     * @brief Computes the sandwich product of a 2D vector and a 2D rotor (R*v*Reverse(R)).
     * @param v The vector to be transformed.
     * @param r The rotor to perform the transformation with.
     * @return The resulting vector after the sandwich product.
     */
    constexpr Vector<T> sandwich(const Vector<T>& v, const Rotor<T>& r) {
        const T s = r.scalar.value;
        const T b = r.bivector.xy;

        // R^2 v, expanded. With e12 v = (y, -x):
        //   R^2 v = (s^2 - b^2) v + 2sb (y, -x).
        // For rotor_from_angle(angle), s^2 - b^2 = cos(angle) and
        // -2sb = sin(angle), so this is the ordinary rotation matrix.
        const T c = s*s - b*b;
        const T d = 2 * s*b;

        return Vector<T>(
            c * v.x + d * v.y,
            c * v.y - d * v.x
        );
    }

    /**
     * @brief Rotates a 2D vector using a 2D rotor.
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
     * @return The bivector, scaled by the rotor's squared norm. The rotor and
     *         the pseudoscalar commute -- the even subalgebra is commutative --
     *         so a unit rotor leaves it exactly unchanged. There is only one
     *         plane in two dimensions, and turning it gives the same plane.
     */
    constexpr Bivector<T> sandwich(const Bivector<T>& b, const Rotor<T>& r) {
        return grade2(r * b * reverse(r));
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
     * @return The rotated bivector, which for a unit rotor is b itself.
     */
    constexpr Bivector<T> rotate(const Bivector<T>& b, const Rotor<T>& r) {
        return sandwich(b, r);
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
} // namespace CliffordCore::Cl2
