#pragma once

/**
 * @file cliffordcore/pga/operations/involutions.hpp
 * @brief The grade involution and Clifford conjugation; reversion is in reverse.hpp.
 *
 * The three involutions, by the sign they give grade k (0 through 4):
 *
 *   reverse    (-1)^(k(k-1)/2)   + + - - +
 *   involute   (-1)^k            + - + - +
 *   conjugate  (-1)^(k(k+1)/2)   + - - + +
 *
 * conjugate == reverse of involute == involute of reverse. None of them
 * depends on the metric, so the signs are the Cl3 ones extended to grade 4.
 * Clifford conjugation is what the general multivector inverse is built
 * on (see inverse.hpp).
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
#include "reverse.hpp"

namespace CliffordCore::PGA
{
    // -----------------------------------------------------------------------
    // Grade involution: negate the odd grades
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Grade involution of a scalar.
     * @param s The scalar.
     * @return The same scalar: grade 0 is even.
     */
    constexpr Scalar<T> involute(const Scalar<T>& s) {
        return s;
    }

    template<typename T>
    /**
     * @brief Grade involution of a plane.
     * @param v The vector.
     * @return The negated vector: grade 1 is odd.
     */
    constexpr Vector<T> involute(const Vector<T>& v) {
        return -v;
    }

    template<typename T>
    /**
     * @brief Grade involution of a line.
     * @param b The bivector.
     * @return The same bivector: grade 2 is even.
     */
    constexpr Bivector<T> involute(const Bivector<T>& b) {
        return b;
    }

    template<typename T>
    /**
     * @brief Grade involution of a point.
     * @param t The trivector.
     * @return The negated trivector: grade 3 is odd.
     */
    constexpr Trivector<T> involute(const Trivector<T>& t) {
        return -t;
    }

    template<typename T>
    /**
     * @brief Grade involution of a quadvector.
     * @param q The quadvector.
     * @return The same quadvector: grade 4 is even.
     */
    constexpr Quadvector<T> involute(const Quadvector<T>& q) {
        return q;
    }

    template<typename T>
    /**
     * @brief Grade involution of a multivector.
     * @param m The multivector.
     * @return The multivector with its odd grades (1 and 3) negated.
     */
    constexpr Multivector<T> involute(const Multivector<T>& m) {
        return Multivector<T>(m.scalar, -m.vector, m.bivector, -m.trivector, m.quadvector);
    }

    template<typename T>
    /**
     * @brief Grade involution of a rotor.
     * @param r The rotor.
     * @return The same rotor: it has only even grades.
     */
    constexpr Rotor<T> involute(const Rotor<T>& r) {
        return r;
    }

    template<typename T>
    /**
     * @brief Grade involution of a translator.
     * @param t The translator.
     * @return The same translator: it has only even grades.
     */
    constexpr Translator<T> involute(const Translator<T>& t) {
        return t;
    }

    template<typename T>
    /**
     * @brief Grade involution of a motor.
     * @param m The motor.
     * @return The same motor: it has only even grades.
     */
    constexpr Motor<T> involute(const Motor<T>& m) {
        return m;
    }

    // -----------------------------------------------------------------------
    // Clifford conjugation: negate grades 1 and 2
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Clifford conjugate of a scalar.
     * @param s The scalar.
     * @return The same scalar.
     */
    constexpr Scalar<T> conjugate(const Scalar<T>& s) {
        return s;
    }

    template<typename T>
    /**
     * @brief Clifford conjugate of a plane.
     * @param v The vector.
     * @return The negated vector.
     */
    constexpr Vector<T> conjugate(const Vector<T>& v) {
        return -v;
    }

    template<typename T>
    /**
     * @brief Clifford conjugate of a line.
     * @param b The bivector.
     * @return The negated bivector.
     */
    constexpr Bivector<T> conjugate(const Bivector<T>& b) {
        return -b;
    }

    template<typename T>
    /**
     * @brief Clifford conjugate of a point.
     * @param t The trivector.
     * @return The same trivector: grade 3 is fixed by conjugation.
     */
    constexpr Trivector<T> conjugate(const Trivector<T>& t) {
        return t;
    }

    template<typename T>
    /**
     * @brief Clifford conjugate of a quadvector.
     * @param q The quadvector.
     * @return The same quadvector: grade 4 is fixed by conjugation.
     */
    constexpr Quadvector<T> conjugate(const Quadvector<T>& q) {
        return q;
    }

    template<typename T>
    /**
     * @brief Clifford conjugate of a multivector.
     * @param m The multivector.
     * @return The multivector with grades 1 and 2 negated. Equal to
     *         reverse(involute(m)), and to involute(reverse(m)).
     */
    constexpr Multivector<T> conjugate(const Multivector<T>& m) {
        return Multivector<T>(m.scalar, -m.vector, -m.bivector, m.trivector, m.quadvector);
    }

    template<typename T>
    /**
     * @brief Clifford conjugate of a rotor.
     * @param r The rotor.
     * @return The rotor with its bivector part negated; for a rotor this is the reverse.
     */
    constexpr Rotor<T> conjugate(const Rotor<T>& r) {
        return reverse(r);
    }

    template<typename T>
    /**
     * @brief Clifford conjugate of a translator.
     * @param t The translator.
     * @return The translator with its ideal part negated; for a translator this is the reverse.
     */
    constexpr Translator<T> conjugate(const Translator<T>& t) {
        return reverse(t);
    }

    template<typename T>
    /**
     * @brief Clifford conjugate of a motor.
     * @param m The motor.
     * @return The motor with its bivector part negated; for a motor this is the reverse.
     */
    constexpr Motor<T> conjugate(const Motor<T>& m) {
        return reverse(m);
    }
} // namespace CliffordCore::PGA
