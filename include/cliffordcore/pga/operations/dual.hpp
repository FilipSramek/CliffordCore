#pragma once

/**
 * @file cliffordcore/pga/operations/dual.hpp
 * @brief Duality as the right complement, and its inverse undual.
 *
 * In Cl3 the dual is multiplication by the pseudoscalar. That cannot work
 * here: e0123 squares to zero, so multiplying by it destroys information
 * and cannot be undone. Cl(3,0,1) uses the complement instead: dual(a) is
 * the blade b with a ^ b = +e0123, extended linearly. It maps grade k to
 * grade 4 - k, swaps the Euclidean and ideal parts, and is what makes the
 * regressive product (join) definable. The signed permutation is:
 *
 *   1     -> +e0123    e01 -> +e23    e012 -> +e3
 *   e0    -> +e123     e02 -> -e13    e013 -> -e2
 *   e1    -> -e023     e03 -> +e12    e023 -> +e1
 *   e2    -> +e013     e12 -> +e03    e123 -> -e0
 *   e3    -> -e012     e13 -> -e02    e0123 -> +1
 *                      e23 -> +e01
 *
 * Applying it twice gives +x on grades 0, 2 and 4 and -x on grades 1 and 3,
 * so undual(x) -- the exact inverse -- is dual(x) on even grades and
 * -dual(x) on odd ones. Every sign here was verified against a model of the
 * algebra before being written down; see docs/pga.md.
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

namespace CliffordCore::PGA
{
    // -----------------------------------------------------------------------
    // dual: the right complement, a ^ dual(a) = +e0123
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Dual of a scalar.
     * @param s The scalar.
     * @return The quadvector s e0123.
     */
    constexpr Quadvector<T> dual(const Scalar<T>& s) {
        return Quadvector<T>(s.value);
    }

    template<typename T>
    /**
     * @brief Dual of a plane: the point it is the complement of.
     * @param v The vector.
     * @return The trivector with e0 -> +e123, e1 -> -e023, e2 -> +e013, e3 -> -e012.
     *         The dual of a Euclidean plane through the origin is an ideal
     *         point (its normal direction); the dual of the ideal plane e0 is
     *         the origin.
     */
    constexpr Trivector<T> dual(const Vector<T>& v) {
        return Trivector<T>(-v.e3, v.e2, -v.e1, v.e0);
    }

    template<typename T>
    /**
     * @brief Dual of a line: the line with its direction and moment swapped.
     * @param b The bivector.
     * @return The bivector with e01 -> +e23, e02 -> -e13, e03 -> +e12,
     *         e12 -> +e03, e13 -> -e02, e23 -> +e01.
     */
    constexpr Bivector<T> dual(const Bivector<T>& b) {
        return Bivector<T>(b.e23, -b.e13, b.e12, b.e03, -b.e02, b.e01);
    }

    template<typename T>
    /**
     * @brief Dual of a point: the plane it is the complement of.
     * @param t The trivector.
     * @return The vector with e012 -> +e3, e013 -> -e2, e023 -> +e1, e123 -> -e0.
     */
    constexpr Vector<T> dual(const Trivector<T>& t) {
        return Vector<T>(-t.e123, t.e023, -t.e013, t.e012);
    }

    template<typename T>
    /**
     * @brief Dual of a quadvector.
     * @param q The quadvector.
     * @return The scalar with the same coefficient.
     */
    constexpr Scalar<T> dual(const Quadvector<T>& q) {
        return Scalar<T>(q.e0123);
    }

    template<typename T>
    /**
     * @brief Dual of a multivector, grade by grade.
     * @param m The multivector.
     * @return The multivector whose grade 4 - k part is the dual of m's grade k part.
     */
    constexpr Multivector<T> dual(const Multivector<T>& m) {
        return Multivector<T>(
            dual(m.quadvector),
            dual(m.trivector),
            dual(m.bivector),
            dual(m.vector),
            dual(m.scalar)
        );
    }

    // -----------------------------------------------------------------------
    // undual: the inverse of dual, so that undual(dual(x)) == x for every x
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Inverse dual of a scalar.
     * @param s The scalar.
     * @return dual(s): the sign is unchanged on even grades.
     */
    constexpr Quadvector<T> undual(const Scalar<T>& s) {
        return dual(s);
    }

    template<typename T>
    /**
     * @brief Inverse dual of a plane.
     * @param v The vector.
     * @return -dual(v): grade 1 flips sign under a double dual.
     */
    constexpr Trivector<T> undual(const Vector<T>& v) {
        return Trivector<T>(v.e3, -v.e2, v.e1, -v.e0);
    }

    template<typename T>
    /**
     * @brief Inverse dual of a line.
     * @param b The bivector.
     * @return dual(b): the sign is unchanged on even grades.
     */
    constexpr Bivector<T> undual(const Bivector<T>& b) {
        return dual(b);
    }

    template<typename T>
    /**
     * @brief Inverse dual of a point.
     * @param t The trivector.
     * @return -dual(t): grade 3 flips sign under a double dual.
     */
    constexpr Vector<T> undual(const Trivector<T>& t) {
        return Vector<T>(t.e123, -t.e023, t.e013, -t.e012);
    }

    template<typename T>
    /**
     * @brief Inverse dual of a quadvector.
     * @param q The quadvector.
     * @return dual(q): the sign is unchanged on even grades.
     */
    constexpr Scalar<T> undual(const Quadvector<T>& q) {
        return dual(q);
    }

    template<typename T>
    /**
     * @brief Inverse dual of a multivector, grade by grade.
     * @param m The multivector.
     * @return The multivector with undual(dual(m)) == m.
     */
    constexpr Multivector<T> undual(const Multivector<T>& m) {
        return Multivector<T>(
            undual(m.quadvector),
            undual(m.trivector),
            undual(m.bivector),
            undual(m.vector),
            undual(m.scalar)
        );
    }
} // namespace CliffordCore::PGA
