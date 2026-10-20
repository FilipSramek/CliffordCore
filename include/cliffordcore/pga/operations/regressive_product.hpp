#pragma once

/**
 * @file cliffordcore/pga/operations/regressive_product.hpp
 * @brief The regressive product (join), and meet as a named spelling of the wedge.
 *
 * In a plane-based algebra the wedge of two objects is their meet: two
 * planes wedge to the line they share, a plane and a line to the point
 * they share. The dual operation -- the object spanned by two others: the
 * line through two points, the plane through a line and a point -- is the
 * regressive product, join(a, b) = undual(dual(a) ^ dual(b)). It is the
 * wedge computed on the complements, which is why Cl(3,0,1) needs the
 * complement-style dual of dual.hpp rather than Cl3's pseudoscalar one.
 *
 * Both are named functions: Cl3 set the precedent of not inventing
 * operators whose reading is not universal, and the & some libraries use
 * for join reads as "and" to most C++ programmers. meet() is the wedge
 * under its geometric name.
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
#include "dual.hpp"
#include "wedge_product.hpp"
#include "mixed_products.hpp"

namespace CliffordCore::PGA
{
    // -----------------------------------------------------------------------
    // Regressive product
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Regressive product of two points: the line through both.
     * @param a The first trivector.
     * @param b The second trivector.
     * @return The resulting bivector. For unit-weight points its Euclidean
     *         part is the direction from a to b and its norm their distance;
     *         it is zero iff the points coincide.
     */
    [[nodiscard]] constexpr Bivector<T> regressive_product(const Trivector<T>& a, const Trivector<T>& b) noexcept {
        return undual(dual(a) ^ dual(b));
    }

    template<typename T>
    /**
     * @brief Regressive product of a point and a line: the plane through both.
     * @param t The trivector.
     * @param b The bivector.
     * @return The resulting vector; zero iff the point lies on the line. Its
     *         norm is the distance from the point to the line, scaled by
     *         their norms.
     */
    [[nodiscard]] constexpr Vector<T> regressive_product(const Trivector<T>& t, const Bivector<T>& b) noexcept {
        return undual(dual(t) ^ dual(b));
    }

    template<typename T>
    /**
     * @brief Regressive product of a line and a point: the plane through both.
     * @param b The bivector.
     * @param t The trivector.
     * @return The resulting vector; zero iff the point lies on the line.
     */
    [[nodiscard]] constexpr Vector<T> regressive_product(const Bivector<T>& b, const Trivector<T>& t) noexcept {
        return undual(dual(b) ^ dual(t));
    }

    template<typename T>
    /**
     * @brief Regressive product of two lines.
     * @param a The first bivector.
     * @param b The second bivector.
     * @return The resulting scalar; zero iff the lines are coplanar.
     */
    [[nodiscard]] constexpr Scalar<T> regressive_product(const Bivector<T>& a, const Bivector<T>& b) noexcept {
        return undual(dual(a) ^ dual(b));
    }

    // -----------------------------------------------------------------------
    // join and meet: the geometric names
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief The line through two points.
     * @param a The first trivector.
     * @param b The second trivector.
     * @return regressive_product(a, b).
     */
    [[nodiscard]] constexpr Bivector<T> join(const Trivector<T>& a, const Trivector<T>& b) noexcept {
        return regressive_product(a, b);
    }

    template<typename T>
    /**
     * @brief The plane through a point and a line.
     * @param t The trivector.
     * @param b The bivector.
     * @return regressive_product(t, b).
     */
    [[nodiscard]] constexpr Vector<T> join(const Trivector<T>& t, const Bivector<T>& b) noexcept {
        return regressive_product(t, b);
    }

    template<typename T>
    /**
     * @brief The plane through a line and a point.
     * @param b The bivector.
     * @param t The trivector.
     * @return regressive_product(b, t).
     */
    [[nodiscard]] constexpr Vector<T> join(const Bivector<T>& b, const Trivector<T>& t) noexcept {
        return regressive_product(b, t);
    }

    template<typename T>
    /**
     * @brief The join of two lines: a scalar that vanishes iff they are coplanar.
     * @param a The first bivector.
     * @param b The second bivector.
     * @return regressive_product(a, b).
     */
    [[nodiscard]] constexpr Scalar<T> join(const Bivector<T>& a, const Bivector<T>& b) noexcept {
        return regressive_product(a, b);
    }

    template<typename T>
    /**
     * @brief The line where two planes meet.
     * @param a The first vector.
     * @param b The second vector.
     * @return a ^ b.
     */
    [[nodiscard]] constexpr Bivector<T> meet(const Vector<T>& a, const Vector<T>& b) noexcept {
        return a ^ b;
    }

    template<typename T>
    /**
     * @brief The point where a plane meets a line.
     * @param v The vector.
     * @param b The bivector.
     * @return v ^ b.
     */
    [[nodiscard]] constexpr Trivector<T> meet(const Vector<T>& v, const Bivector<T>& b) noexcept {
        return v ^ b;
    }

    template<typename T>
    /**
     * @brief The point where a line meets a plane.
     * @param b The bivector.
     * @param v The vector.
     * @return b ^ v.
     */
    [[nodiscard]] constexpr Trivector<T> meet(const Bivector<T>& b, const Vector<T>& v) noexcept {
        return b ^ v;
    }
} // namespace CliffordCore::PGA
