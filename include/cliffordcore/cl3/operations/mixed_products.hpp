#pragma once

/**
 * @file cliffordcore/cl3/operations/mixed_products.hpp
 * @brief Geometric and wedge products between differing grades.
 *
 * Products between differing grades: each promotes both operands to
 * Multivector, defers to the Cayley table, and returns the most general
 * type the product can produce. operator* is defined for every one of the
 * 36 operand pairs.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include "addition.hpp"
#include "geometric_product.hpp"

// Geometric products between differing grades.
//
// Every one of these promotes both operands to a multivector, defers to the
// general product in geometric_product.hpp, and returns a Multivector. That
// keeps them correct by construction -- there is one multiplication table in the
// library and these all go through it -- and matches the rule used elsewhere:
// operator* returns the most general type the product can produce.
//
// Not defined here: Vector * Vector, Multivector * Multivector and
// Rotor * Rotor have direct implementations in geometric_product.hpp, and the
// Scalar pairings are members of the individual types.
//
// Trivector * Trivector is not here either: two pseudoscalars always multiply
// to a pure scalar, so Trivector defines it directly and returns Scalar<T>.
// The tests assert it agrees with the general product's scalar part.

namespace CliffordCore::Cl3
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
            return Multivector<T>(r.scalar, Vector<T>(), r.bivector, Trivector<T>());
        }
    } // namespace detail

    template<typename T>
    /**
     * @brief Geometric product of a vector and a bivector.
     * @param v The vector operand.
     * @param b The bivector operand.
     * @return The resulting multivector, carrying vector and trivector parts.
     */
    constexpr Multivector<T> operator*(const Vector<T>& v, const Bivector<T>& b) {
        return geometric_product(detail::promote(v), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a bivector and a vector.
     * @param b The bivector operand.
     * @param v The vector operand.
     * @return The resulting multivector, carrying vector and trivector parts.
     */
    constexpr Multivector<T> operator*(const Bivector<T>& b, const Vector<T>& v) {
        return geometric_product(detail::promote(b), detail::promote(v));
    }

    template<typename T>
    /**
     * @brief Geometric product of two bivectors.
     * @param a The left bivector.
     * @param b The right bivector.
     * @return The resulting multivector, carrying scalar and bivector parts.
     */
    constexpr Multivector<T> operator*(const Bivector<T>& a, const Bivector<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a vector and a trivector.
     * @param v The vector operand.
     * @param t The trivector operand.
     * @return The resulting multivector, carrying a bivector part.
     */
    constexpr Multivector<T> operator*(const Vector<T>& v, const Trivector<T>& t) {
        return geometric_product(detail::promote(v), detail::promote(t));
    }

    template<typename T>
    /**
     * @brief Geometric product of a trivector and a vector.
     * @param t The trivector operand.
     * @param v The vector operand.
     * @return The resulting multivector, carrying a bivector part.
     */
    constexpr Multivector<T> operator*(const Trivector<T>& t, const Vector<T>& v) {
        return geometric_product(detail::promote(t), detail::promote(v));
    }

    template<typename T>
    /**
     * @brief Geometric product of a bivector and a trivector.
     * @param b The bivector operand.
     * @param t The trivector operand.
     * @return The resulting multivector, carrying a vector part.
     */
    constexpr Multivector<T> operator*(const Bivector<T>& b, const Trivector<T>& t) {
        return geometric_product(detail::promote(b), detail::promote(t));
    }

    template<typename T>
    /**
     * @brief Geometric product of a trivector and a bivector.
     * @param t The trivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector, carrying a vector part.
     */
    constexpr Multivector<T> operator*(const Trivector<T>& t, const Bivector<T>& b) {
        return geometric_product(detail::promote(t), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a rotor and a vector.
     * @param r The rotor operand.
     * @param v The vector operand.
     * @return The resulting multivector. This is the piece needed to spell a
     *         rotation as r * v * reverse(r).
     */
    constexpr Multivector<T> operator*(const Rotor<T>& r, const Vector<T>& v) {
        return geometric_product(detail::promote(r), detail::promote(v));
    }

    template<typename T>
    /**
     * @brief Geometric product of a vector and a rotor.
     * @param v The vector operand.
     * @param r The rotor operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Vector<T>& v, const Rotor<T>& r) {
        return geometric_product(detail::promote(v), detail::promote(r));
    }

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
     * @brief Geometric product of a multivector and a trivector.
     * @param m The multivector operand.
     * @param t The trivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Multivector<T>& m, const Trivector<T>& t) {
        return geometric_product(m, detail::promote(t));
    }

    template<typename T>
    /**
     * @brief Geometric product of a trivector and a multivector.
     * @param t The trivector operand.
     * @param m The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Trivector<T>& t, const Multivector<T>& m) {
        return geometric_product(detail::promote(t), m);
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

    template<typename T>
    /**
     * @brief Geometric product of a rotor and a bivector.
     * @param r The rotor operand.
     * @param b The bivector operand.
     * @return The resulting multivector. Needed to sandwich a bivector.
     */
    constexpr Multivector<T> operator*(const Rotor<T>& r, const Bivector<T>& b) {
        return geometric_product(detail::promote(r), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a bivector and a rotor.
     * @param b The bivector operand.
     * @param r The rotor operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Bivector<T>& b, const Rotor<T>& r) {
        return geometric_product(detail::promote(b), detail::promote(r));
    }

    template<typename T>
    /**
     * @brief Geometric product of a rotor and a trivector.
     * @param r The rotor operand.
     * @param t The trivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Rotor<T>& r, const Trivector<T>& t) {
        return geometric_product(detail::promote(r), detail::promote(t));
    }

    template<typename T>
    /**
     * @brief Geometric product of a trivector and a rotor.
     * @param t The trivector operand.
     * @param r The rotor operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Trivector<T>& t, const Rotor<T>& r) {
        return geometric_product(detail::promote(t), detail::promote(r));
    }

    // -----------------------------------------------------------------------
    // Wedge beyond grade 1 with grade 1
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Wedge product of a vector and a bivector, giving the pseudoscalar.
     * @param v The vector operand.
     * @param b The bivector operand.
     * @return The resulting trivector, the grade 3 part of their product. This
     *         is the only geometric route to a trivector in this library.
     */
    constexpr Trivector<T> operator^(const Vector<T>& v, const Bivector<T>& b) {
        return geometric_product(detail::promote(v), detail::promote(b)).trivector;
    }

    template<typename T>
    /**
     * @brief Wedge product of a bivector and a vector, giving the pseudoscalar.
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
     * @brief Wedge product of a vector and a bivector.
     * @param v The vector operand.
     * @param b The bivector operand.
     * @return The resulting trivector.
     */
    constexpr Trivector<T> wedge_product(const Vector<T>& v, const Bivector<T>& b) {
        return v ^ b;
    }

    template<typename T>
    /**
     * @brief Wedge product of a bivector and a vector.
     * @param b The bivector operand.
     * @param v The vector operand.
     * @return The resulting trivector.
     */
    constexpr Trivector<T> wedge_product(const Bivector<T>& b, const Vector<T>& v) {
        return b ^ v;
    }
} // namespace CliffordCore::Cl3
