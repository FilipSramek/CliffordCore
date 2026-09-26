#pragma once

/**
 * @file cliffordcore/cl2/operations/reverse.hpp
 * @brief The reverse involution, which flips grade 2.
 *
 * Reverse scales grade k by (-1)^(k(k-1)/2): grades 0 and 1 keep their
 * sign, grade 2 is negated. It is the involution rotors are undone with,
 * since inverse(r) is reverse(r) for a unit rotor. Read as a complex number,
 * a rotor's reverse is its complex conjugate.
 *
 * @author Filip Sramek
 * @version 0.3.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"

namespace CliffordCore::Cl2
{
    template<typename T>
    /**
     * @brief Computes the geometric reverse of a scalar.
     * @param s The scalar for which to compute the reverse.
     * @return The scalar unchanged. Grade 0 keeps its sign under reverse.
     */
    constexpr Scalar<T> reverse(const Scalar<T>& s) {
        return s;
    }

    template<typename T>
    /**
     * @brief Computes the geometric reverse of a 2D vector.
     * @param v The vector for which to compute the reverse.
     * @return The reverse of the vector v. (reverse(v) = v)
     */
    constexpr Vector<T> reverse(const Vector<T>& v) {
        return Vector<T>(v.x, v.y);
    }

    template<typename T>
    /**
     * @brief Computes the geometric reverse of a bivector in 2D space.
     * @param b The bivector for which to compute the reverse.
     * @return The reverse of the bivector b. (reverse(b) = -b)
     */
    constexpr Bivector<T> reverse(const Bivector<T>& b) {
        return Bivector<T>(-b.xy);
    }

    template<typename T>
    /**
     * @brief Computes the geometric reverse of a multivector in 2D space.
     * @param m The multivector for which to compute the reverse.
     * @return The reverse of the multivector m. (reverse(m) = scalar(m) + vector(m) - bivector(m))
     */
    constexpr Multivector<T> reverse(const Multivector<T>& m) {
        return Multivector<T>(
            m.scalar,
            m.vector,
            Bivector<T>(-m.bivector.xy)
        );
    }

    template<typename T>
    /**
     * @brief Computes the geometric reverse of a rotor in 2D space.
     * @param r The rotor for which to compute the reverse.
     * @return The reverse of the rotor r. (reverse(r) = scalar(r) - bivector(r))
     */
    constexpr Rotor<T> reverse(const Rotor<T>& r) {
        return Rotor<T>(r.scalar, Bivector<T>(-r.bivector.xy));
    }
} // namespace CliffordCore::Cl2
