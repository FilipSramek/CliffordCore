#pragma once

/**
 * @file cliffordcore/pga/operations/mixed_products.hpp
 * @brief Geometric and wedge products between differing types.
 *
 * Products between differing types: each promotes both operands to
 * Multivector, defers to the Cayley table, and returns the most general
 * type the product can produce. Together with the type headers and
 * geometric_product.hpp this makes operator* defined for every one of the
 * 81 operand pairs.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include "addition.hpp"
#include "geometric_product.hpp"

// Geometric products between differing types.
//
// Every one of these promotes both operands to a multivector, defers to the
// general product in geometric_product.hpp, and returns the most general type
// the product can produce: a Motor when both operands are even (the even
// subalgebra is closed under the product, so nothing is lost) and a
// Multivector otherwise. That keeps them correct by construction -- there is
// one multiplication table in the library and these all go through it.
//
// Not defined here: Vector * Vector, Multivector * Multivector, Rotor * Rotor,
// Translator * Translator and Motor * Motor have direct implementations in
// geometric_product.hpp; the Scalar pairings are members of the individual
// types; and Quadvector * Quadvector is a member of Quadvector (it is always
// zero).
//
// Unlike Cl3, Trivector * Trivector IS here and returns a Multivector: two
// points share exactly two basis indices, so their product has grades 0 and 2.

namespace CliffordCore::PGA
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
            return Multivector<T>(r);
        }

        template<typename T>
        /**
         * @brief Widens a translator into a multivector.
         * @param t The translator to widen.
         * @return The equivalent multivector.
         */
        constexpr Multivector<T> promote(const Translator<T>& t) {
            return Multivector<T>(t);
        }

        template<typename T>
        /**
         * @brief Widens a motor into a multivector.
         * @param m The motor to widen.
         * @return The equivalent multivector.
         */
        constexpr Multivector<T> promote(const Motor<T>& m) {
            return Multivector<T>(m);
        }
    } // namespace detail

    // -----------------------------------------------------------------------
    // Pairs that can leave the even subalgebra: every one returns Multivector
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Geometric product of a vector (plane) and a bivector (line).
     * @param a The vector operand.
     * @param b The bivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Vector<T>& a, const Bivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a vector (plane) and a trivector (point).
     * @param a The vector operand.
     * @param b The trivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Vector<T>& a, const Trivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a vector (plane) and a quadvector.
     * @param a The vector operand.
     * @param b The quadvector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Vector<T>& a, const Quadvector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a vector (plane) and a multivector.
     * @param a The vector operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Vector<T>& a, const Multivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a vector (plane) and a rotor.
     * @param a The vector operand.
     * @param b The rotor operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Vector<T>& a, const Rotor<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a vector (plane) and a translator.
     * @param a The vector operand.
     * @param b The translator operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Vector<T>& a, const Translator<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a vector (plane) and a motor.
     * @param a The vector operand.
     * @param b The motor operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Vector<T>& a, const Motor<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a bivector (line) and a vector (plane).
     * @param a The bivector operand.
     * @param b The vector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Bivector<T>& a, const Vector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a bivector (line) and a bivector (line).
     * @param a The bivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector, carrying grades 0, 2 and 4.
     */
    constexpr Multivector<T> operator*(const Bivector<T>& a, const Bivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a bivector (line) and a trivector (point).
     * @param a The bivector operand.
     * @param b The trivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Bivector<T>& a, const Trivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a bivector (line) and a quadvector.
     * @param a The bivector operand.
     * @param b The quadvector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Bivector<T>& a, const Quadvector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a bivector (line) and a multivector.
     * @param a The bivector operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Bivector<T>& a, const Multivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a bivector (line) and a rotor.
     * @param a The bivector operand.
     * @param b The rotor operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Bivector<T>& a, const Rotor<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a bivector (line) and a translator.
     * @param a The bivector operand.
     * @param b The translator operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Bivector<T>& a, const Translator<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a bivector (line) and a motor.
     * @param a The bivector operand.
     * @param b The motor operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Bivector<T>& a, const Motor<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a trivector (point) and a vector (plane).
     * @param a The trivector operand.
     * @param b The vector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Trivector<T>& a, const Vector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a trivector (point) and a bivector (line).
     * @param a The trivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Trivector<T>& a, const Bivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a trivector (point) and a trivector (point).
     * @param a The trivector operand.
     * @param b The trivector operand.
     * @return The resulting multivector, carrying grades 0 and 2 -- two points multiply to a scalar plus an ideal line, not to a scalar as in Cl3.
     */
    constexpr Multivector<T> operator*(const Trivector<T>& a, const Trivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a trivector (point) and a quadvector.
     * @param a The trivector operand.
     * @param b The quadvector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Trivector<T>& a, const Quadvector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a trivector (point) and a multivector.
     * @param a The trivector operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Trivector<T>& a, const Multivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a trivector (point) and a rotor.
     * @param a The trivector operand.
     * @param b The rotor operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Trivector<T>& a, const Rotor<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a trivector (point) and a translator.
     * @param a The trivector operand.
     * @param b The translator operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Trivector<T>& a, const Translator<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a trivector (point) and a motor.
     * @param a The trivector operand.
     * @param b The motor operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Trivector<T>& a, const Motor<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a quadvector and a vector (plane).
     * @param a The quadvector operand.
     * @param b The vector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Quadvector<T>& a, const Vector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a quadvector and a bivector (line).
     * @param a The quadvector operand.
     * @param b The bivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Quadvector<T>& a, const Bivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a quadvector and a trivector (point).
     * @param a The quadvector operand.
     * @param b The trivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Quadvector<T>& a, const Trivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a quadvector and a multivector.
     * @param a The quadvector operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Quadvector<T>& a, const Multivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a quadvector and a rotor.
     * @param a The quadvector operand.
     * @param b The rotor operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Quadvector<T>& a, const Rotor<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a quadvector and a translator.
     * @param a The quadvector operand.
     * @param b The translator operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Quadvector<T>& a, const Translator<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a quadvector and a motor.
     * @param a The quadvector operand.
     * @param b The motor operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Quadvector<T>& a, const Motor<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a multivector and a vector (plane).
     * @param a The multivector operand.
     * @param b The vector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Multivector<T>& a, const Vector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a multivector and a bivector (line).
     * @param a The multivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Multivector<T>& a, const Bivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a multivector and a trivector (point).
     * @param a The multivector operand.
     * @param b The trivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Multivector<T>& a, const Trivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a multivector and a quadvector.
     * @param a The multivector operand.
     * @param b The quadvector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Multivector<T>& a, const Quadvector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a multivector and a rotor.
     * @param a The multivector operand.
     * @param b The rotor operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Multivector<T>& a, const Rotor<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a multivector and a translator.
     * @param a The multivector operand.
     * @param b The translator operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Multivector<T>& a, const Translator<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a multivector and a motor.
     * @param a The multivector operand.
     * @param b The motor operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Multivector<T>& a, const Motor<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a rotor and a vector (plane).
     * @param a The rotor operand.
     * @param b The vector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Rotor<T>& a, const Vector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a rotor and a bivector (line).
     * @param a The rotor operand.
     * @param b The bivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Rotor<T>& a, const Bivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a rotor and a trivector (point).
     * @param a The rotor operand.
     * @param b The trivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Rotor<T>& a, const Trivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a rotor and a quadvector.
     * @param a The rotor operand.
     * @param b The quadvector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Rotor<T>& a, const Quadvector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a rotor and a multivector.
     * @param a The rotor operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Rotor<T>& a, const Multivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a translator and a vector (plane).
     * @param a The translator operand.
     * @param b The vector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Translator<T>& a, const Vector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a translator and a bivector (line).
     * @param a The translator operand.
     * @param b The bivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Translator<T>& a, const Bivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a translator and a trivector (point).
     * @param a The translator operand.
     * @param b The trivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Translator<T>& a, const Trivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a translator and a quadvector.
     * @param a The translator operand.
     * @param b The quadvector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Translator<T>& a, const Quadvector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a translator and a multivector.
     * @param a The translator operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Translator<T>& a, const Multivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a motor and a vector (plane).
     * @param a The motor operand.
     * @param b The vector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Motor<T>& a, const Vector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a motor and a bivector (line).
     * @param a The motor operand.
     * @param b The bivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Motor<T>& a, const Bivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a motor and a trivector (point).
     * @param a The motor operand.
     * @param b The trivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Motor<T>& a, const Trivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a motor and a quadvector.
     * @param a The motor operand.
     * @param b The quadvector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Motor<T>& a, const Quadvector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a motor and a multivector.
     * @param a The motor operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Motor<T>& a, const Multivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    // -----------------------------------------------------------------------
    // Pairs of even types: the product stays even, so these return Motor
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Geometric product of a rotor and a translator.
     * @param a The rotor operand.
     * @param b The translator operand.
     * @return The resulting motor, equivalent to applying b then a. The even
     *         subalgebra is closed, so narrowing the general product loses nothing.
     */
    constexpr Motor<T> operator*(const Rotor<T>& a, const Translator<T>& b) {
        return Motor<T>(geometric_product(detail::promote(a), detail::promote(b)));
    }

    template<typename T>
    /**
     * @brief Geometric product of a rotor and a motor.
     * @param a The rotor operand.
     * @param b The motor operand.
     * @return The resulting motor, equivalent to applying b then a. The even
     *         subalgebra is closed, so narrowing the general product loses nothing.
     */
    constexpr Motor<T> operator*(const Rotor<T>& a, const Motor<T>& b) {
        return Motor<T>(geometric_product(detail::promote(a), detail::promote(b)));
    }

    template<typename T>
    /**
     * @brief Geometric product of a translator and a rotor.
     * @param a The translator operand.
     * @param b The rotor operand.
     * @return The resulting motor, equivalent to applying b then a. The even
     *         subalgebra is closed, so narrowing the general product loses nothing.
     */
    constexpr Motor<T> operator*(const Translator<T>& a, const Rotor<T>& b) {
        return Motor<T>(geometric_product(detail::promote(a), detail::promote(b)));
    }

    template<typename T>
    /**
     * @brief Geometric product of a translator and a motor.
     * @param a The translator operand.
     * @param b The motor operand.
     * @return The resulting motor, equivalent to applying b then a. The even
     *         subalgebra is closed, so narrowing the general product loses nothing.
     */
    constexpr Motor<T> operator*(const Translator<T>& a, const Motor<T>& b) {
        return Motor<T>(geometric_product(detail::promote(a), detail::promote(b)));
    }

    template<typename T>
    /**
     * @brief Geometric product of a motor and a rotor.
     * @param a The motor operand.
     * @param b The rotor operand.
     * @return The resulting motor, equivalent to applying b then a. The even
     *         subalgebra is closed, so narrowing the general product loses nothing.
     */
    constexpr Motor<T> operator*(const Motor<T>& a, const Rotor<T>& b) {
        return Motor<T>(geometric_product(detail::promote(a), detail::promote(b)));
    }

    template<typename T>
    /**
     * @brief Geometric product of a motor and a translator.
     * @param a The motor operand.
     * @param b The translator operand.
     * @return The resulting motor, equivalent to applying b then a. The even
     *         subalgebra is closed, so narrowing the general product loses nothing.
     */
    constexpr Motor<T> operator*(const Motor<T>& a, const Translator<T>& b) {
        return Motor<T>(geometric_product(detail::promote(a), detail::promote(b)));
    }

    template<typename T>
    /**
     * @brief Geometric product of two points, packed as a motor.
     * @param a The first trivector.
     * @param b The second trivector.
     * @return The resulting motor: the composition of the point reflections
     *         in b and then a, which for unit-weight points is minus the
     *         translator by twice the vector from b to a. Same components as
     *         a * b, whose odd grades are zero.
     */
    constexpr Motor<T> motor_product(const Trivector<T>& a, const Trivector<T>& b) {
        return Motor<T>(geometric_product(detail::promote(a), detail::promote(b)));
    }

    // -----------------------------------------------------------------------
    // Wedge (meet) beyond two planes
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Wedge product of a plane and a line: the point where they meet.
     * @param v The vector operand.
     * @param b The bivector operand.
     * @return The resulting trivector, the grade 3 part of their product.
     */
    constexpr Trivector<T> operator^(const Vector<T>& v, const Bivector<T>& b) {
        return geometric_product(detail::promote(v), detail::promote(b)).trivector;
    }

    template<typename T>
    /**
     * @brief Wedge product of a line and a plane: the point where they meet.
     * @param b The bivector operand.
     * @param v The vector operand.
     * @return The resulting trivector. The wedge is graded-commutative here, so
     *         this equals v ^ b.
     */
    constexpr Trivector<T> operator^(const Bivector<T>& b, const Vector<T>& v) {
        return geometric_product(detail::promote(b), detail::promote(v)).trivector;
    }

    template<typename T>
    /**
     * @brief Wedge product of a plane and a point: the incidence test.
     * @param v The vector operand.
     * @param t The trivector operand.
     * @return The resulting quadvector. For plane(a, b, c, d) and point(x, y, z)
     *         its coefficient is a x + b y + c z + d, so it is zero iff the
     *         point lies on the plane.
     */
    constexpr Quadvector<T> operator^(const Vector<T>& v, const Trivector<T>& t) {
        return geometric_product(detail::promote(v), detail::promote(t)).quadvector;
    }

    template<typename T>
    /**
     * @brief Wedge product of a point and a plane: the incidence test.
     * @param t The trivector operand.
     * @param v The vector operand.
     * @return The resulting quadvector, equal to -(v ^ t).
     */
    constexpr Quadvector<T> operator^(const Trivector<T>& t, const Vector<T>& v) {
        return geometric_product(detail::promote(t), detail::promote(v)).quadvector;
    }

    template<typename T>
    /**
     * @brief Wedge product of two lines.
     * @param a The first bivector.
     * @param b The second bivector.
     * @return The resulting quadvector. Zero iff the two lines are coplanar
     *         (they meet or are parallel); non-zero for skew lines.
     */
    constexpr Quadvector<T> operator^(const Bivector<T>& a, const Bivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b)).quadvector;
    }

    template<typename T>
    /**
     * @brief Wedge product of a plane and a line.
     * @param v The vector operand.
     * @param b The bivector operand.
     * @return The resulting trivector.
     */
    constexpr Trivector<T> wedge_product(const Vector<T>& v, const Bivector<T>& b) {
        return v ^ b;
    }

    template<typename T>
    /**
     * @brief Wedge product of a line and a plane.
     * @param b The bivector operand.
     * @param v The vector operand.
     * @return The resulting trivector.
     */
    constexpr Trivector<T> wedge_product(const Bivector<T>& b, const Vector<T>& v) {
        return b ^ v;
    }

    template<typename T>
    /**
     * @brief Wedge product of a plane and a point.
     * @param v The vector operand.
     * @param t The trivector operand.
     * @return The resulting quadvector.
     */
    constexpr Quadvector<T> wedge_product(const Vector<T>& v, const Trivector<T>& t) {
        return v ^ t;
    }

    template<typename T>
    /**
     * @brief Wedge product of a point and a plane.
     * @param t The trivector operand.
     * @param v The vector operand.
     * @return The resulting quadvector.
     */
    constexpr Quadvector<T> wedge_product(const Trivector<T>& t, const Vector<T>& v) {
        return t ^ v;
    }

    template<typename T>
    /**
     * @brief Wedge product of two lines.
     * @param a The first bivector.
     * @param b The second bivector.
     * @return The resulting quadvector.
     */
    constexpr Quadvector<T> wedge_product(const Bivector<T>& a, const Bivector<T>& b) {
        return a ^ b;
    }
} // namespace CliffordCore::PGA
