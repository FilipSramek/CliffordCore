#pragma once

/**
 * @file cliffordcore/cl2/operations/contraction.hpp
 * @brief Left and right contractions, and the scalar product.
 *
 * operator| is the left contraction. For two vectors that is exactly the
 * dot product, so nothing familiar changes; the contraction is simply what
 * it generalises to. Cl(2,0) has only three grades, so the family is small:
 * vector into bivector, and bivector with bivector.
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
// negative the result is zero, so `Bivector _| Vector` is identically zero and
// is not provided -- writing it is almost always a mistake for `Vector _| B`.
//
// Cl(2,0) makes these unusually simple. A vector times a bivector is ALREADY a
// pure vector, and a bivector times a bivector is already a pure scalar, so the
// grade filter selects the whole product in every case below. The contractions
// are still spelled out, and still named, so that code written against Cl(3,0)
// reads the same here.
//
// Scalar contractions are just scaling (`s _| A == s * A`, `A _| s == 0` for
// grade above 0), so they are omitted here; use the existing operator*.

namespace CliffordCore::Cl2
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
     *         Defined wherever a*b has a grade 0 part: not for a vector with a
     *         bivector or a rotor, whose product is a pure vector.
     */
    constexpr auto scalar_product(const A& a, const B& b) -> decltype(grade0(a * b)) {
        return grade0(a * b);
    }

    // -----------------------------------------------------------------------
    // Left contraction
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Left contraction of two vectors.
     * @param a The left vector.
     * @param b The right vector.
     * @return The resulting scalar. For equal grades this is the dot product.
     */
    constexpr Scalar<T> left_contraction(const Vector<T>& a, const Vector<T>& b) {
        return a | b;
    }

    template<typename T>
    /**
     * @brief Left contraction of a vector into a bivector.
     * @param v The vector operand.
     * @param b The bivector operand.
     * @return The resulting vector: v turned a quarter turn counter-clockwise
     *         and scaled by b. This is the product the projection and
     *         rejection formulas are written with.
     */
    constexpr Vector<T> left_contraction(const Vector<T>& v, const Bivector<T>& b) {
        return v * b;
    }

    template<typename T>
    /**
     * @brief Left contraction of two bivectors.
     * @param a The left bivector.
     * @param b The right bivector.
     * @return The resulting scalar, -a.xy * b.xy.
     */
    constexpr Scalar<T> left_contraction(const Bivector<T>& a, const Bivector<T>& b) {
        return a * b;
    }

    // -----------------------------------------------------------------------
    // operator| -- the left contraction. Vector | Vector lives in
    // dot_product.hpp and is unchanged; these extend it across grades.
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Left contraction of a vector into a bivector.
     * @param v The vector operand.
     * @param b The bivector operand.
     * @return The resulting vector.
     */
    constexpr Vector<T> operator|(const Vector<T>& v, const Bivector<T>& b) {
        return left_contraction(v, b);
    }

    template<typename T>
    /**
     * @brief Left contraction of two bivectors.
     * @param a The left bivector.
     * @param b The right bivector.
     * @return The resulting scalar.
     */
    constexpr Scalar<T> operator|(const Bivector<T>& a, const Bivector<T>& b) {
        return left_contraction(a, b);
    }

    // -----------------------------------------------------------------------
    // Right contraction: A |_ B = reverse(reverse(B) _| reverse(A)).
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Right contraction of two vectors.
     * @param a The left vector.
     * @param b The right vector.
     * @return The resulting scalar, equal to the left contraction here.
     */
    constexpr Scalar<T> right_contraction(const Vector<T>& a, const Vector<T>& b) {
        return a | b;
    }

    template<typename T>
    /**
     * @brief Right contraction of a bivector by a vector.
     * @param b The bivector operand.
     * @param v The vector operand.
     * @return The resulting vector. Because e12 anticommutes with vectors, this
     *         is the negation of left_contraction(v, b).
     */
    constexpr Vector<T> right_contraction(const Bivector<T>& b, const Vector<T>& v) {
        return b * v;
    }

    template<typename T>
    /**
     * @brief Right contraction of two bivectors.
     * @param a The left bivector.
     * @param b The right bivector.
     * @return The resulting scalar, equal to the left contraction here.
     */
    constexpr Scalar<T> right_contraction(const Bivector<T>& a, const Bivector<T>& b) {
        return a * b;
    }
} // namespace CliffordCore::Cl2
