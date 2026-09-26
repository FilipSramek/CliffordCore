#pragma once

/**
 * @file cliffordcore/cl2/operations/inverse.hpp
 * @brief Multiplicative inverses for every type.
 *
 * Inverses carry the minus signs the squares imply: v^-1 = v/|v|^2, but
 * B^-1 = -B/|B|^2. A Multivector inverts through Clifford conjugation, and
 * in Cl(2,0) that is a single division: m * conjugate(m) is a pure scalar.
 *
 * @author Filip Sramek
 * @version 0.3.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <type_traits>
#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"
#include "norm.hpp"
#include "reverse.hpp"
#include "involutions.hpp"

namespace CliffordCore::Cl2
{
    template<typename T>
    /**
     * @brief Computes the inverse of a scalar.
     * @param s The scalar for which to compute the inverse.
     * @return The inverse scalar to the input scalar.
     */
    constexpr Scalar<T> inverse(const Scalar<T>& s) {
        return Scalar<T>(T(1) / s.value);
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a 2D vector.
     * @param v The vector for which to compute the inverse.
     * @return The inverse vector to the input vector.
     */
    constexpr Vector<T> inverse(const Vector<T>& v) {
        return v/squared_norm(v);
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a bivector in 2D space.
     * @param b The bivector for which to compute the inverse.
     * @return The inverse bivector to the input bivector. The pseudoscalar
     *         squares to -1, so the inverse carries a minus sign: b^-1 = -b/|b|^2.
     */
    constexpr Bivector<T> inverse(const Bivector<T>& b) {
        return -b/squared_norm(b);
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a multivector in 2D space.
     * @param m The multivector for which to compute the inverse.
     * @return The inverse multivector, satisfying m * inverse(m) = 1.
     *
     * Dividing by the squared norm only works for a single grade. In general,
     * with m = a + x e1 + y e2 + b e12, the product m * conjugate(m) is
     *
     *     a^2 - x^2 - y^2 + b^2
     *
     * and nothing else. The cross terms v*B + B*v cancel because the
     * pseudoscalar anticommutes with vectors in two dimensions, so -- unlike
     * Cl(3,0), where a pseudoscalar part survives and has to be inverted like
     * a complex number -- this is a pure scalar, and
     * inverse(m) = conjugate(m) / (a^2 - x^2 - y^2 + b^2).
     *
     * Mind the minus signs: that denominator is not a sum of squares, so a
     * non-zero multivector can still be singular. 1 + e1 is the standard
     * example, since (1 + e1)(1 - e1) = 1 - e1^2 = 0. When the denominator is
     * zero this returns a zero multivector rather than dividing by zero.
     */
    constexpr Multivector<T> inverse(const Multivector<T>& m) {
        const T a = m.scalar.value;
        const T x = m.vector.x;
        const T y = m.vector.y;
        const T b = m.bivector.xy;

        const T denominator = a * a - x * x - y * y + b * b;

        if (denominator == T(0)) {
            return Multivector<T>();
        }

        return conjugate(m) / Scalar<T>(denominator);
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a rotor in 2D space.
     * @param r The rotor for which to compute the inverse.
     * @return The inverse rotor to the input rotor: r^-1 = reverse(r)/|r|^2,
     *         since r * reverse(r) = |r|^2. For a unit rotor this is just the
     *         reverse, which is why rotations undo by reversing. Read as a
     *         complex number, this is the familiar conj(z)/|z|^2.
     */
    constexpr Rotor<T> inverse(const Rotor<T>& r) {
        return reverse(r)/squared_norm(r);
    }
} // namespace CliffordCore::Cl2
