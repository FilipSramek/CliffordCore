#pragma once

/**
 * @file cliffordcore/cl2/operations/involutions.hpp
 * @brief The grade involution and Clifford conjugation.
 *
 * The grade involution negates the odd grade, 1; Clifford conjugation
 * negates grades 1 and 2. Together with reverse they form the set of three,
 * and conjugate is exactly reverse composed with involute.
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
#include "reverse.hpp"

// Cl(2,0) has three sign-flipping involutions, and they come as a set. Each
// scales grade k by a sign that depends only on k:
//
//              grade:            0    1    2
//   reverse    (-1)^(k(k-1)/2)   +    +    -
//   involute   (-1)^k            +    -    +
//   conjugate  (-1)^(k(k+1)/2)   +    -    -
//
// So conjugate == reverse composed with involute, in either order. The test
// suite asserts that identity rather than trusting the table.
//
// reverse() lives in reverse.hpp, which this header includes, so including
// this one gives you all three.

namespace CliffordCore::Cl2
{
    template<typename T>
    /**
     * @brief Applies the grade involution to a scalar.
     * @param s The scalar to involve.
     * @return The scalar unchanged; grade 0 keeps its sign.
     */
    constexpr Scalar<T> involute(const Scalar<T>& s) {
        return s;
    }

    template<typename T>
    /**
     * @brief Applies the grade involution to a vector.
     * @param v The vector to involve.
     * @return The negated vector; grade 1 flips sign.
     */
    constexpr Vector<T> involute(const Vector<T>& v) {
        return Vector<T>(-v.x, -v.y);
    }

    template<typename T>
    /**
     * @brief Applies the grade involution to a bivector.
     * @param b The bivector to involve.
     * @return The bivector unchanged; grade 2 keeps its sign.
     */
    constexpr Bivector<T> involute(const Bivector<T>& b) {
        return b;
    }

    template<typename T>
    /**
     * @brief Applies the grade involution to a multivector.
     * @param m The multivector to involve.
     * @return The multivector with its odd grade, 1, negated.
     */
    constexpr Multivector<T> involute(const Multivector<T>& m) {
        return Multivector<T>(
            m.scalar,
            Vector<T>(-m.vector.x, -m.vector.y),
            m.bivector
        );
    }

    template<typename T>
    /**
     * @brief Applies the grade involution to a rotor.
     * @param r The rotor to involve.
     * @return The rotor unchanged. A rotor holds only even grades, which is
     *         exactly what makes rotors closed under the sandwich product.
     */
    constexpr Rotor<T> involute(const Rotor<T>& r) {
        return r;
    }

    // -----------------------------------------------------------------------
    // Clifford conjugation: reverse composed with the grade involution.
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Applies Clifford conjugation to a scalar.
     * @param s The scalar to conjugate.
     * @return The scalar unchanged.
     */
    constexpr Scalar<T> conjugate(const Scalar<T>& s) {
        return s;
    }

    template<typename T>
    /**
     * @brief Applies Clifford conjugation to a vector.
     * @param v The vector to conjugate.
     * @return The negated vector.
     */
    constexpr Vector<T> conjugate(const Vector<T>& v) {
        return Vector<T>(-v.x, -v.y);
    }

    template<typename T>
    /**
     * @brief Applies Clifford conjugation to a bivector.
     * @param b The bivector to conjugate.
     * @return The negated bivector.
     */
    constexpr Bivector<T> conjugate(const Bivector<T>& b) {
        return Bivector<T>(-b.xy);
    }

    template<typename T>
    /**
     * @brief Applies Clifford conjugation to a multivector.
     * @param m The multivector to conjugate.
     * @return The multivector with grades 1 and 2 negated.
     *
     * This is the involution that makes the general inverse work. In Cl(2,0)
     * m multiplied by its conjugate collapses all the way to a pure scalar --
     * no pseudoscalar part survives, unlike Cl(3,0) -- so inverting m is a
     * single division. See inverse.hpp.
     */
    constexpr Multivector<T> conjugate(const Multivector<T>& m) {
        return Multivector<T>(
            m.scalar,
            Vector<T>(-m.vector.x, -m.vector.y),
            Bivector<T>(-m.bivector.xy)
        );
    }

    template<typename T>
    /**
     * @brief Applies Clifford conjugation to a rotor.
     * @param r The rotor to conjugate.
     * @return The rotor with its bivector part negated. For a rotor this
     *         coincides with reverse(), since a rotor has no odd grades.
     */
    constexpr Rotor<T> conjugate(const Rotor<T>& r) {
        return Rotor<T>(r.scalar, Bivector<T>(-r.bivector.xy));
    }
} // namespace CliffordCore::Cl2
