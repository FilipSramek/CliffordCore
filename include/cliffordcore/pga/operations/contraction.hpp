#pragma once

/**
 * @file cliffordcore/pga/operations/contraction.hpp
 * @brief Left and right contractions, and the scalar product.
 *
 * operator| is the left contraction, exactly as in Cl3. Cl3 chose it over
 * the symmetric "fat dot" because it is the inner product that survives a
 * degenerate metric -- and this is that metric. Every contraction here is
 * a grade projection of the general product, so e0^2 = 0 is accounted for
 * automatically: contracting anything into an ideal element, or an ideal
 * element into anything, silently loses the e0 terms as it must.
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
#include "addition.hpp"
#include "dot_product.hpp"
#include "geometric_product.hpp"
#include "grade.hpp"
#include "mixed_products.hpp"

// The inner-product family.
//
// For blades A of grade r and B of grade s, the LEFT contraction is the part of
// their geometric product with grade s - r:
//
//     A _| B  =  <A B>_(s-r)
//
// and the RIGHT contraction takes grade r - s instead. Where the grade would go
// negative the result is zero, so those pairs are not provided. Bivectors in
// Cl(3,0,1) need not be blades, but grade-filtering the general product is
// still the definition of the contraction on a sum of blades, so it stays
// correct; the tests check every overload against the product.
//
// The geometric reading: a plane contracted into a point is the line through
// the point perpendicular to the plane, a plane into a line is the plane
// through the line perpendicular to the first, and so on -- which is why the
// projections in geometry.hpp are written with these.
//
// Scalar contractions are just scaling and are omitted; use operator*.

namespace CliffordCore::PGA
{
    // -----------------------------------------------------------------------
    // Scalar product: the grade 0 part of any product.
    // -----------------------------------------------------------------------

    template<typename A, typename B>
    /**
     * @brief Computes the scalar product of two elements.
     * @param a The left operand.
     * @param b The right operand.
     * @return The grade 0 part of a*b, the one piece every product shares.
     */
    [[nodiscard]] constexpr auto scalar_product(const A& a, const B& b) noexcept -> decltype(grade0(a * b)) {
        return grade0(a * b);
    }

    // -----------------------------------------------------------------------
    // Left contraction
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Left contraction of two planes.
     * @param a The left vector.
     * @param b The right vector.
     * @return The resulting scalar. For equal grades this is the dot product.
     */
    [[nodiscard]] constexpr Scalar<T> left_contraction(const Vector<T>& a, const Vector<T>& b) noexcept {
        return a | b;
    }

    template<typename T>
    /**
     * @brief Left contraction of a plane into a line.
     * @param v The vector.
     * @param b The bivector.
     * @return The grade 1 part of v*b: the plane through the line perpendicular to v.
     */
    [[nodiscard]] constexpr Vector<T> left_contraction(const Vector<T>& v, const Bivector<T>& b) noexcept {
        return grade1(v * b);
    }

    template<typename T>
    /**
     * @brief Left contraction of a plane into a point.
     * @param v The vector.
     * @param t The trivector.
     * @return The grade 2 part of v*t: the line through the point perpendicular to v.
     */
    [[nodiscard]] constexpr Bivector<T> left_contraction(const Vector<T>& v, const Trivector<T>& t) noexcept {
        return grade2(v * t);
    }

    template<typename T>
    /**
     * @brief Left contraction of a plane into the pseudoscalar.
     * @param v The vector.
     * @param q The quadvector.
     * @return The grade 3 part of v*q: an ideal point, the plane's normal direction.
     */
    [[nodiscard]] constexpr Trivector<T> left_contraction(const Vector<T>& v, const Quadvector<T>& q) noexcept {
        return grade3(v * q);
    }

    template<typename T>
    /**
     * @brief Left contraction of two lines.
     * @param a The left bivector.
     * @param b The right bivector.
     * @return The grade 0 part of a*b: minus the dot product of their directions.
     */
    [[nodiscard]] constexpr Scalar<T> left_contraction(const Bivector<T>& a, const Bivector<T>& b) noexcept {
        return grade0(a * b);
    }

    template<typename T>
    /**
     * @brief Left contraction of a line into a point.
     * @param b The bivector.
     * @param t The trivector.
     * @return The grade 1 part of b*t: the plane through the point perpendicular to b.
     */
    [[nodiscard]] constexpr Vector<T> left_contraction(const Bivector<T>& b, const Trivector<T>& t) noexcept {
        return grade1(b * t);
    }

    template<typename T>
    /**
     * @brief Left contraction of a line into the pseudoscalar.
     * @param b The bivector.
     * @param q The quadvector.
     * @return The grade 2 part of b*q: an ideal line.
     */
    [[nodiscard]] constexpr Bivector<T> left_contraction(const Bivector<T>& b, const Quadvector<T>& q) noexcept {
        return grade2(b * q);
    }

    template<typename T>
    /**
     * @brief Left contraction of two points.
     * @param a The left trivector.
     * @param b The right trivector.
     * @return The grade 0 part of a*b: minus the product of their weights.
     */
    [[nodiscard]] constexpr Scalar<T> left_contraction(const Trivector<T>& a, const Trivector<T>& b) noexcept {
        return grade0(a * b);
    }

    template<typename T>
    /**
     * @brief Left contraction of a point into the pseudoscalar.
     * @param t The trivector.
     * @param q The quadvector.
     * @return The grade 1 part of t*q: an ideal plane.
     */
    [[nodiscard]] constexpr Vector<T> left_contraction(const Trivector<T>& t, const Quadvector<T>& q) noexcept {
        return grade1(t * q);
    }

    template<typename T>
    /**
     * @brief Left contraction of two quadvectors.
     * @param a The left quadvector.
     * @param b The right quadvector.
     * @return Zero, always: e0123 squares to zero.
     */
    [[nodiscard]] constexpr Scalar<T> left_contraction(const Quadvector<T>& a, const Quadvector<T>& b) noexcept {
        return a * b;
    }

    // operator| is the left contraction. Vector | Vector lives in dot_product.hpp.

    template<typename T>
    /**
     * @brief Left contraction of a plane into a line.
     * @param v The vector.
     * @param b The bivector.
     * @return left_contraction(v, b).
     */
    [[nodiscard]] constexpr Vector<T> operator|(const Vector<T>& v, const Bivector<T>& b) noexcept {
        return left_contraction(v, b);
    }

    template<typename T>
    /**
     * @brief Left contraction of a plane into a point.
     * @param v The vector.
     * @param t The trivector.
     * @return left_contraction(v, t).
     */
    [[nodiscard]] constexpr Bivector<T> operator|(const Vector<T>& v, const Trivector<T>& t) noexcept {
        return left_contraction(v, t);
    }

    template<typename T>
    /**
     * @brief Left contraction of a plane into the pseudoscalar.
     * @param v The vector.
     * @param q The quadvector.
     * @return left_contraction(v, q).
     */
    [[nodiscard]] constexpr Trivector<T> operator|(const Vector<T>& v, const Quadvector<T>& q) noexcept {
        return left_contraction(v, q);
    }

    template<typename T>
    /**
     * @brief Left contraction of two lines.
     * @param a The left bivector.
     * @param b The right bivector.
     * @return left_contraction(a, b).
     */
    [[nodiscard]] constexpr Scalar<T> operator|(const Bivector<T>& a, const Bivector<T>& b) noexcept {
        return left_contraction(a, b);
    }

    template<typename T>
    /**
     * @brief Left contraction of a line into a point.
     * @param b The bivector.
     * @param t The trivector.
     * @return left_contraction(b, t).
     */
    [[nodiscard]] constexpr Vector<T> operator|(const Bivector<T>& b, const Trivector<T>& t) noexcept {
        return left_contraction(b, t);
    }

    template<typename T>
    /**
     * @brief Left contraction of a line into the pseudoscalar.
     * @param b The bivector.
     * @param q The quadvector.
     * @return left_contraction(b, q).
     */
    [[nodiscard]] constexpr Bivector<T> operator|(const Bivector<T>& b, const Quadvector<T>& q) noexcept {
        return left_contraction(b, q);
    }

    template<typename T>
    /**
     * @brief Left contraction of two points.
     * @param a The left trivector.
     * @param b The right trivector.
     * @return left_contraction(a, b).
     */
    [[nodiscard]] constexpr Scalar<T> operator|(const Trivector<T>& a, const Trivector<T>& b) noexcept {
        return left_contraction(a, b);
    }

    template<typename T>
    /**
     * @brief Left contraction of a point into the pseudoscalar.
     * @param t The trivector.
     * @param q The quadvector.
     * @return left_contraction(t, q).
     */
    [[nodiscard]] constexpr Vector<T> operator|(const Trivector<T>& t, const Quadvector<T>& q) noexcept {
        return left_contraction(t, q);
    }

    template<typename T>
    /**
     * @brief Left contraction of two quadvectors.
     * @param a The left quadvector.
     * @param b The right quadvector.
     * @return left_contraction(a, b), which is zero.
     */
    [[nodiscard]] constexpr Scalar<T> operator|(const Quadvector<T>& a, const Quadvector<T>& b) noexcept {
        return left_contraction(a, b);
    }

    // -----------------------------------------------------------------------
    // Right contraction: the mirror image, grade r - s of A B.
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Right contraction of two planes.
     * @param a The left vector.
     * @param b The right vector.
     * @return The resulting scalar, equal to the left contraction for equal grades.
     */
    [[nodiscard]] constexpr Scalar<T> right_contraction(const Vector<T>& a, const Vector<T>& b) noexcept {
        return a | b;
    }

    template<typename T>
    /**
     * @brief Right contraction of a line by a plane.
     * @param b The bivector.
     * @param v The vector.
     * @return The grade 1 part of b*v.
     */
    [[nodiscard]] constexpr Vector<T> right_contraction(const Bivector<T>& b, const Vector<T>& v) noexcept {
        return grade1(b * v);
    }

    template<typename T>
    /**
     * @brief Right contraction of a point by a plane.
     * @param t The trivector.
     * @param v The vector.
     * @return The grade 2 part of t*v.
     */
    [[nodiscard]] constexpr Bivector<T> right_contraction(const Trivector<T>& t, const Vector<T>& v) noexcept {
        return grade2(t * v);
    }

    template<typename T>
    /**
     * @brief Right contraction of the pseudoscalar by a plane.
     * @param q The quadvector.
     * @param v The vector.
     * @return The grade 3 part of q*v.
     */
    [[nodiscard]] constexpr Trivector<T> right_contraction(const Quadvector<T>& q, const Vector<T>& v) noexcept {
        return grade3(q * v);
    }

    template<typename T>
    /**
     * @brief Right contraction of two lines.
     * @param a The left bivector.
     * @param b The right bivector.
     * @return The grade 0 part of a*b.
     */
    [[nodiscard]] constexpr Scalar<T> right_contraction(const Bivector<T>& a, const Bivector<T>& b) noexcept {
        return grade0(a * b);
    }

    template<typename T>
    /**
     * @brief Right contraction of a point by a line.
     * @param t The trivector.
     * @param b The bivector.
     * @return The grade 1 part of t*b.
     */
    [[nodiscard]] constexpr Vector<T> right_contraction(const Trivector<T>& t, const Bivector<T>& b) noexcept {
        return grade1(t * b);
    }

    template<typename T>
    /**
     * @brief Right contraction of the pseudoscalar by a line.
     * @param q The quadvector.
     * @param b The bivector.
     * @return The grade 2 part of q*b.
     */
    [[nodiscard]] constexpr Bivector<T> right_contraction(const Quadvector<T>& q, const Bivector<T>& b) noexcept {
        return grade2(q * b);
    }

    template<typename T>
    /**
     * @brief Right contraction of two points.
     * @param a The left trivector.
     * @param b The right trivector.
     * @return The grade 0 part of a*b.
     */
    [[nodiscard]] constexpr Scalar<T> right_contraction(const Trivector<T>& a, const Trivector<T>& b) noexcept {
        return grade0(a * b);
    }

    template<typename T>
    /**
     * @brief Right contraction of the pseudoscalar by a point.
     * @param q The quadvector.
     * @param t The trivector.
     * @return The grade 1 part of q*t.
     */
    [[nodiscard]] constexpr Vector<T> right_contraction(const Quadvector<T>& q, const Trivector<T>& t) noexcept {
        return grade1(q * t);
    }

    template<typename T>
    /**
     * @brief Right contraction of two quadvectors.
     * @param a The left quadvector.
     * @param b The right quadvector.
     * @return Zero, always.
     */
    [[nodiscard]] constexpr Scalar<T> right_contraction(const Quadvector<T>& a, const Quadvector<T>& b) noexcept {
        return a * b;
    }
} // namespace CliffordCore::PGA
