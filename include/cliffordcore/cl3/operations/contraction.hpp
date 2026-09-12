#pragma once

/**
 * @file cliffordcore/cl3/operations/contraction.hpp
 * @brief Left and right contractions, and the scalar product.
 *
 * operator| is the left contraction. For two vectors that is exactly the
 * dot product, so nothing familiar changes; the contraction is simply what
 * it generalises to. It is chosen over the symmetric fat dot because it
 * survives a degenerate metric, which will matter when Cl(3,0,1) lands.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"
#include "addition.hpp"
#include "geometric_product.hpp"
#include "grade.hpp"
#include "mixed_products.hpp"
#include "reverse.hpp"

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
// Every bivector in 3D is a blade, so grade-filtering the general product is a
// correct implementation here; the tests check it against that product.
//
// operator| is the left contraction. For two vectors that is exactly the dot
// product it has always been, so nothing about the existing behaviour changes --
// this only extends it to the mixed-grade cases. The left contraction is chosen
// over the symmetric "fat dot" deliberately: it is the one that keeps working in
// a degenerate metric, which matters for Cl(3,0,1) later, and it is the product
// the projection formulas are written with.
//
// Scalar contractions are just scaling (`s _| A == s * A`, `A _| s == 0` for
// grade above 0), so they are omitted here; use the existing operator*.

namespace CliffordCore::Cl3
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
     * @return The resulting vector: the part of b perpendicular to v, within b.
     *         This is the product the projection onto a plane is written with.
     */
    constexpr Vector<T> left_contraction(const Vector<T>& v, const Bivector<T>& b) {
        return grade1(v * b);
    }

    template<typename T>
    /**
     * @brief Left contraction of a vector into a trivector.
     * @param v The vector operand.
     * @param t The trivector operand.
     * @return The resulting bivector.
     */
    constexpr Bivector<T> left_contraction(const Vector<T>& v, const Trivector<T>& t) {
        return grade2(v * t);
    }

    template<typename T>
    /**
     * @brief Left contraction of two bivectors.
     * @param a The left bivector.
     * @param b The right bivector.
     * @return The resulting scalar.
     */
    constexpr Scalar<T> left_contraction(const Bivector<T>& a, const Bivector<T>& b) {
        return grade0(a * b);
    }

    template<typename T>
    /**
     * @brief Left contraction of a bivector into a trivector.
     * @param b The bivector operand.
     * @param t The trivector operand.
     * @return The resulting vector.
     */
    constexpr Vector<T> left_contraction(const Bivector<T>& b, const Trivector<T>& t) {
        return grade1(b * t);
    }

    template<typename T>
    /**
     * @brief Left contraction of two trivectors.
     * @param a The left trivector.
     * @param b The right trivector.
     * @return The resulting scalar.
     */
    constexpr Scalar<T> left_contraction(const Trivector<T>& a, const Trivector<T>& b) {
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
     * @brief Left contraction of a vector into a trivector.
     * @param v The vector operand.
     * @param t The trivector operand.
     * @return The resulting bivector.
     */
    constexpr Bivector<T> operator|(const Vector<T>& v, const Trivector<T>& t) {
        return left_contraction(v, t);
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

    template<typename T>
    /**
     * @brief Left contraction of a bivector into a trivector.
     * @param b The bivector operand.
     * @param t The trivector operand.
     * @return The resulting vector.
     */
    constexpr Vector<T> operator|(const Bivector<T>& b, const Trivector<T>& t) {
        return left_contraction(b, t);
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
     * @return The resulting vector.
     */
    constexpr Vector<T> right_contraction(const Bivector<T>& b, const Vector<T>& v) {
        return grade1(b * v);
    }

    template<typename T>
    /**
     * @brief Right contraction of a trivector by a vector.
     * @param t The trivector operand.
     * @param v The vector operand.
     * @return The resulting bivector.
     */
    constexpr Bivector<T> right_contraction(const Trivector<T>& t, const Vector<T>& v) {
        return grade2(t * v);
    }

    template<typename T>
    /**
     * @brief Right contraction of two bivectors.
     * @param a The left bivector.
     * @param b The right bivector.
     * @return The resulting scalar.
     */
    constexpr Scalar<T> right_contraction(const Bivector<T>& a, const Bivector<T>& b) {
        return grade0(a * b);
    }

    template<typename T>
    /**
     * @brief Right contraction of a trivector by a bivector.
     * @param t The trivector operand.
     * @param b The bivector operand.
     * @return The resulting vector.
     */
    constexpr Vector<T> right_contraction(const Trivector<T>& t, const Bivector<T>& b) {
        return grade1(t * b);
    }

    template<typename T>
    /**
     * @brief Right contraction of two trivectors.
     * @param a The left trivector.
     * @param b The right trivector.
     * @return The resulting scalar.
     */
    constexpr Scalar<T> right_contraction(const Trivector<T>& a, const Trivector<T>& b) {
        return a * b;
    }
} // namespace CliffordCore::Cl3
