#pragma once

/**
 * @file cliffordcore/cl2/operations/mixed_products.hpp
 * @brief Geometric products between differing grades.
 *
 * Every product here promotes both operands, defers to the general product
 * in geometric_product.hpp, and returns the most general type the product
 * can produce. In Cl(2,0) that is often narrower than a Multivector: a
 * vector times a bivector is a pure vector, and a rotor times a vector is a
 * pure vector too.
 *
 * @author Filip Sramek
 * @version 0.3.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include "addition.hpp"
#include "geometric_product.hpp"

// Geometric products between differing grades.
//
// Every one of these promotes both operands to a multivector and defers to the
// general product in geometric_product.hpp. That keeps them correct by
// construction -- there is one multiplication table in the library and these
// all go through it -- and matches the rule used elsewhere: operator* returns
// the most general type the product can produce.
//
// In Cl(2,0) that rule gives a TIGHTER table than in Cl(3,0). The odd part of
// the algebra is grade 1 alone and the even part is grades 0 and 2, so:
//
//   odd  * even  -> odd   : Vector * Bivector, Vector * Rotor    -> Vector
//   even * odd   -> odd   : Bivector * Vector, Rotor * Vector    -> Vector
//   even * even  -> even  : Bivector * Rotor, Rotor * Bivector   -> Rotor
//
// Those six return the narrow type, taken as one grade (or the even pair) of
// the general product -- the other grades are zero by construction, not by
// approximation. Everything involving a Multivector still returns a
// Multivector.
//
// Not defined here: Vector * Vector, Multivector * Multivector and
// Rotor * Rotor have direct implementations in geometric_product.hpp,
// Bivector * Bivector is a member of Bivector (it returns a Scalar, since
// e12^2 = -1), and the Scalar pairings are members of the individual types.
//
// There is no wedge beyond grade 1 with grade 1: Vector ^ Bivector would be
// grade 3, which does not exist here. It is identically zero, and left
// undefined rather than returning a type that pretends otherwise.

namespace CliffordCore::Cl2
{
    namespace detail
    {
        template<typename T>
        /**
         * @brief Widens a rotor into a multivector.
         * @param r The rotor to widen.
         * @return The equivalent multivector.
         */
        constexpr Multivector<T> promote(const Rotor<T>& r) {
            return Multivector<T>(r.scalar, Vector<T>(), r.bivector);
        }
    } // namespace detail

    // -----------------------------------------------------------------------
    // Odd times even, and even times odd: always a pure vector
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Geometric product of a vector and a bivector.
     * @param v The vector operand.
     * @param b The bivector operand.
     * @return The resulting vector. With b = beta e12 this is v turned a
     *         quarter turn counter-clockwise and scaled by beta: e1 e12 = e2,
     *         e2 e12 = -e1. It is dual(v) when beta = 1.
     */
    constexpr Vector<T> operator*(const Vector<T>& v, const Bivector<T>& b) {
        return geometric_product(detail::promote(v), detail::promote(b)).vector;
    }

    template<typename T>
    /**
     * @brief Geometric product of a bivector and a vector.
     * @param b The bivector operand.
     * @param v The vector operand.
     * @return The resulting vector. The pseudoscalar anticommutes with vectors,
     *         so this is the negation of v * b.
     */
    constexpr Vector<T> operator*(const Bivector<T>& b, const Vector<T>& v) {
        return geometric_product(detail::promote(b), detail::promote(v)).vector;
    }

    template<typename T>
    /**
     * @brief Geometric product of a vector and a rotor.
     * @param v The vector operand.
     * @param r The rotor operand.
     * @return The resulting vector.
     */
    constexpr Vector<T> operator*(const Vector<T>& v, const Rotor<T>& r) {
        return geometric_product(detail::promote(v), detail::promote(r)).vector;
    }

    template<typename T>
    /**
     * @brief Geometric product of a rotor and a vector.
     * @param r The rotor operand.
     * @param v The vector operand.
     * @return The resulting vector. In two dimensions this one-sided product
     *         already rotates: r * r * v is the same as the sandwich
     *         r * v * reverse(r). See sandwich.hpp.
     */
    constexpr Vector<T> operator*(const Rotor<T>& r, const Vector<T>& v) {
        return geometric_product(detail::promote(r), detail::promote(v)).vector;
    }

    // -----------------------------------------------------------------------
    // Even times even: always a rotor
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Geometric product of a bivector and a rotor.
     * @param b The bivector operand.
     * @param r The rotor operand.
     * @return The resulting rotor. The even subalgebra is closed, so the vector
     *         part of the general product is always zero here.
     */
    constexpr Rotor<T> operator*(const Bivector<T>& b, const Rotor<T>& r) {
        const Multivector<T> product = geometric_product(detail::promote(b), detail::promote(r));
        return Rotor<T>(product.scalar, product.bivector);
    }

    template<typename T>
    /**
     * @brief Geometric product of a rotor and a bivector.
     * @param r The rotor operand.
     * @param b The bivector operand.
     * @return The resulting rotor. Needed to sandwich a bivector.
     */
    constexpr Rotor<T> operator*(const Rotor<T>& r, const Bivector<T>& b) {
        const Multivector<T> product = geometric_product(detail::promote(r), detail::promote(b));
        return Rotor<T>(product.scalar, product.bivector);
    }

    // -----------------------------------------------------------------------
    // Anything with a multivector: a multivector
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Geometric product of a multivector and a vector.
     * @param m The multivector operand.
     * @param v The vector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Multivector<T>& m, const Vector<T>& v) {
        return geometric_product(m, detail::promote(v));
    }

    template<typename T>
    /**
     * @brief Geometric product of a vector and a multivector.
     * @param v The vector operand.
     * @param m The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Vector<T>& v, const Multivector<T>& m) {
        return geometric_product(detail::promote(v), m);
    }

    template<typename T>
    /**
     * @brief Geometric product of a multivector and a bivector.
     * @param m The multivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Multivector<T>& m, const Bivector<T>& b) {
        return geometric_product(m, detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a bivector and a multivector.
     * @param b The bivector operand.
     * @param m The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Bivector<T>& b, const Multivector<T>& m) {
        return geometric_product(detail::promote(b), m);
    }

    template<typename T>
    /**
     * @brief Geometric product of a multivector and a rotor.
     * @param m The multivector operand.
     * @param r The rotor operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Multivector<T>& m, const Rotor<T>& r) {
        return geometric_product(m, detail::promote(r));
    }

    template<typename T>
    /**
     * @brief Geometric product of a rotor and a multivector.
     * @param r The rotor operand.
     * @param m The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Rotor<T>& r, const Multivector<T>& m) {
        return geometric_product(detail::promote(r), m);
    }
} // namespace CliffordCore::Cl2
