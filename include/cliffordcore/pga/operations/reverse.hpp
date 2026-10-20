#pragma once

/**
 * @file cliffordcore/pga/operations/reverse.hpp
 * @brief Reversion: reversing the order of the basis vectors in each blade.
 *
 * Reversion flips grade k by (-1)^(k(k-1)/2): grades 0, 1 and 4 keep their
 * sign, grades 2 and 3 are negated. It does not depend on the metric, so
 * the pattern is the same as Cl3's with grade 4 added. For a unit motor
 * the reverse is the inverse, which is why motions undo by reversing.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../quadvector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"
#include "../translator.hpp"
#include "../motor.hpp"

namespace CliffordCore::PGA
{
    template<typename T>
    /**
     * @brief Reverse of a scalar.
     * @param s The scalar.
     * @return The same scalar: grade 0 is unchanged by reversion.
     */
    [[nodiscard]] constexpr Scalar<T> reverse(const Scalar<T>& s) noexcept {
        return s;
    }

    template<typename T>
    /**
     * @brief Reverse of a plane.
     * @param v The vector.
     * @return The same vector: grade 1 is unchanged by reversion.
     */
    [[nodiscard]] constexpr Vector<T> reverse(const Vector<T>& v) noexcept {
        return v;
    }

    template<typename T>
    /**
     * @brief Reverse of a line.
     * @param b The bivector.
     * @return The negated bivector: grade 2 changes sign under reversion.
     */
    [[nodiscard]] constexpr Bivector<T> reverse(const Bivector<T>& b) noexcept {
        return -b;
    }

    template<typename T>
    /**
     * @brief Reverse of a point.
     * @param t The trivector.
     * @return The negated trivector: grade 3 changes sign under reversion.
     */
    [[nodiscard]] constexpr Trivector<T> reverse(const Trivector<T>& t) noexcept {
        return -t;
    }

    template<typename T>
    /**
     * @brief Reverse of a quadvector.
     * @param q The quadvector.
     * @return The same quadvector: grade 4 is unchanged by reversion.
     */
    [[nodiscard]] constexpr Quadvector<T> reverse(const Quadvector<T>& q) noexcept {
        return q;
    }

    template<typename T>
    /**
     * @brief Reverse of a multivector.
     * @param m The multivector.
     * @return The multivector with its grade 2 and grade 3 parts negated.
     */
    [[nodiscard]] constexpr Multivector<T> reverse(const Multivector<T>& m) noexcept {
        return Multivector<T>(m.scalar, m.vector, -m.bivector, -m.trivector, m.quadvector);
    }

    template<typename T>
    /**
     * @brief Reverse of a rotor.
     * @param r The rotor.
     * @return The rotor with its bivector part negated: the inverse rotation for a unit rotor.
     */
    [[nodiscard]] constexpr Rotor<T> reverse(const Rotor<T>& r) noexcept {
        return Rotor<T>(r.scalar, -r.e12, -r.e13, -r.e23);
    }

    template<typename T>
    /**
     * @brief Reverse of a translator.
     * @param t The translator.
     * @return The translator with its ideal part negated: the opposite translation for a normalised translator.
     */
    [[nodiscard]] constexpr Translator<T> reverse(const Translator<T>& t) noexcept {
        return Translator<T>(t.scalar, -t.e01, -t.e02, -t.e03);
    }

    template<typename T>
    /**
     * @brief Reverse of a motor.
     * @param m The motor.
     * @return The motor with its bivector part negated and e0123 kept: the inverse motion for a unit motor.
     */
    [[nodiscard]] constexpr Motor<T> reverse(const Motor<T>& m) noexcept {
        return Motor<T>(m.scalar, -m.bivector, m.quadvector);
    }
} // namespace CliffordCore::PGA
