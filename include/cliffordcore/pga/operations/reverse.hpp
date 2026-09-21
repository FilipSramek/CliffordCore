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
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
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
    constexpr Scalar<T> reverse(const Scalar<T>& s) {
        return s;
    }

    template<typename T>
    /**
     * @brief Reverse of a plane.
     * @param v The vector.
     * @return The same vector: grade 1 is unchanged by reversion.
     */
    constexpr Vector<T> reverse(const Vector<T>& v) {
        return v;
    }

    template<typename T>
    /**
     * @brief Reverse of a line.
     * @param b The bivector.
     * @return The negated bivector: grade 2 changes sign under reversion.
     */
    constexpr Bivector<T> reverse(const Bivector<T>& b) {
        return -b;
    }

    template<typename T>
    /**
     * @brief Reverse of a point.
     * @param t The trivector.
     * @return The negated trivector: grade 3 changes sign under reversion.
     */
    constexpr Trivector<T> reverse(const Trivector<T>& t) {
        return -t;
    }

    template<typename T>
    /**
     * @brief Reverse of a quadvector.
     * @param q The quadvector.
     * @return The same quadvector: grade 4 is unchanged by reversion.
     */
    constexpr Quadvector<T> reverse(const Quadvector<T>& q) {
        return q;
    }

    template<typename T>
    /**
     * @brief Reverse of a multivector.
     * @param m The multivector.
     * @return The multivector with its grade 2 and grade 3 parts negated.
     */
    constexpr Multivector<T> reverse(const Multivector<T>& m) {
        return Multivector<T>(m.scalar, m.vector, -m.bivector, -m.trivector, m.quadvector);
    }

    template<typename T>
    /**
     * @brief Reverse of a rotor.
     * @param r The rotor.
     * @return The rotor with its bivector part negated: the inverse rotation for a unit rotor.
     */
    constexpr Rotor<T> reverse(const Rotor<T>& r) {
        return Rotor<T>(r.scalar, -r.e12, -r.e13, -r.e23);
    }

    template<typename T>
    /**
     * @brief Reverse of a translator.
     * @param t The translator.
     * @return The translator with its ideal part negated: the opposite translation for a normalised translator.
     */
    constexpr Translator<T> reverse(const Translator<T>& t) {
        return Translator<T>(t.scalar, -t.e01, -t.e02, -t.e03);
    }

    template<typename T>
    /**
     * @brief Reverse of a motor.
     * @param m The motor.
     * @return The motor with its bivector part negated and e0123 kept: the inverse motion for a unit motor.
     */
    constexpr Motor<T> reverse(const Motor<T>& m) {
        return Motor<T>(m.scalar, -m.bivector, m.quadvector);
    }
} // namespace CliffordCore::PGA
