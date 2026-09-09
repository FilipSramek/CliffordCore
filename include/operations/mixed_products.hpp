#pragma once

/**
 * @file mixed_products.hpp
 * @brief Geometric and wedge products between differing grades.
 */

#include "addition.hpp"
#include "geometric_product.hpp"

// Geometric products between differing grades.
//
// Every one of these promotes both operands to a multivector, defers to the
// general product in geometric_product.hpp, and returns a Multivector3. That
// keeps them correct by construction -- there is one multiplication table in the
// library and these all go through it -- and matches the rule used elsewhere:
// operator* returns the most general type the product can produce.
//
// Not defined here: Vector3 * Vector3, Multivector3 * Multivector3 and
// Rotor3 * Rotor3 have direct implementations in geometric_product.hpp, and the
// Scalar pairings are members of the individual types.
//
// Trivector3 * Trivector3 is not here either: two pseudoscalars always multiply
// to a pure scalar, so Trivector3 defines it directly and returns Scalar<T>.
// The tests assert it agrees with the general product's scalar part.

namespace CliffordCore
{
    namespace detail
    {
        template<typename T>
        /**
         * @brief Widens a rotor into a multivector.
         * @param r The rotor to widen.
         * @return The equivalent multivector.
         */
        constexpr Multivector3<T> promote(const Rotor3<T>& r) {
            return Multivector3<T>(r.scalar, Vector3<T>(), r.bivector, Trivector3<T>());
        }
    } // namespace detail

    template<typename T>
    /**
     * @brief Geometric product of a vector and a bivector.
     * @param v The vector operand.
     * @param b The bivector operand.
     * @return The resulting multivector, carrying vector and trivector parts.
     */
    constexpr Multivector3<T> operator*(const Vector3<T>& v, const Bivector3<T>& b) {
        return geometric_product(detail::promote(v), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a bivector and a vector.
     * @param b The bivector operand.
     * @param v The vector operand.
     * @return The resulting multivector, carrying vector and trivector parts.
     */
    constexpr Multivector3<T> operator*(const Bivector3<T>& b, const Vector3<T>& v) {
        return geometric_product(detail::promote(b), detail::promote(v));
    }

    template<typename T>
    /**
     * @brief Geometric product of two bivectors.
     * @param a The left bivector.
     * @param b The right bivector.
     * @return The resulting multivector, carrying scalar and bivector parts.
     */
    constexpr Multivector3<T> operator*(const Bivector3<T>& a, const Bivector3<T>& b) {
        return geometric_product(detail::promote(a), detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a vector and a trivector.
     * @param v The vector operand.
     * @param t The trivector operand.
     * @return The resulting multivector, carrying a bivector part.
     */
    constexpr Multivector3<T> operator*(const Vector3<T>& v, const Trivector3<T>& t) {
        return geometric_product(detail::promote(v), detail::promote(t));
    }

    template<typename T>
    /**
     * @brief Geometric product of a trivector and a vector.
     * @param t The trivector operand.
     * @param v The vector operand.
     * @return The resulting multivector, carrying a bivector part.
     */
    constexpr Multivector3<T> operator*(const Trivector3<T>& t, const Vector3<T>& v) {
        return geometric_product(detail::promote(t), detail::promote(v));
    }

    template<typename T>
    /**
     * @brief Geometric product of a bivector and a trivector.
     * @param b The bivector operand.
     * @param t The trivector operand.
     * @return The resulting multivector, carrying a vector part.
     */
    constexpr Multivector3<T> operator*(const Bivector3<T>& b, const Trivector3<T>& t) {
        return geometric_product(detail::promote(b), detail::promote(t));
    }

    template<typename T>
    /**
     * @brief Geometric product of a trivector and a bivector.
     * @param t The trivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector, carrying a vector part.
     */
    constexpr Multivector3<T> operator*(const Trivector3<T>& t, const Bivector3<T>& b) {
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
    constexpr Multivector3<T> operator*(const Rotor3<T>& r, const Vector3<T>& v) {
        return geometric_product(detail::promote(r), detail::promote(v));
    }

    template<typename T>
    /**
     * @brief Geometric product of a vector and a rotor.
     * @param v The vector operand.
     * @param r The rotor operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator*(const Vector3<T>& v, const Rotor3<T>& r) {
        return geometric_product(detail::promote(v), detail::promote(r));
    }

    template<typename T>
    /**
     * @brief Geometric product of a multivector and a vector.
     * @param m The multivector operand.
     * @param v The vector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator*(const Multivector3<T>& m, const Vector3<T>& v) {
        return geometric_product(m, detail::promote(v));
    }

    template<typename T>
    /**
     * @brief Geometric product of a vector and a multivector.
     * @param v The vector operand.
     * @param m The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator*(const Vector3<T>& v, const Multivector3<T>& m) {
        return geometric_product(detail::promote(v), m);
    }

    template<typename T>
    /**
     * @brief Geometric product of a multivector and a bivector.
     * @param m The multivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator*(const Multivector3<T>& m, const Bivector3<T>& b) {
        return geometric_product(m, detail::promote(b));
    }

    template<typename T>
    /**
     * @brief Geometric product of a bivector and a multivector.
     * @param b The bivector operand.
     * @param m The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator*(const Bivector3<T>& b, const Multivector3<T>& m) {
        return geometric_product(detail::promote(b), m);
    }

    template<typename T>
    /**
     * @brief Geometric product of a multivector and a trivector.
     * @param m The multivector operand.
     * @param t The trivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator*(const Multivector3<T>& m, const Trivector3<T>& t) {
        return geometric_product(m, detail::promote(t));
    }

    template<typename T>
    /**
     * @brief Geometric product of a trivector and a multivector.
     * @param t The trivector operand.
     * @param m The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator*(const Trivector3<T>& t, const Multivector3<T>& m) {
        return geometric_product(detail::promote(t), m);
    }

    template<typename T>
    /**
     * @brief Geometric product of a multivector and a rotor.
     * @param m The multivector operand.
     * @param r The rotor operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator*(const Multivector3<T>& m, const Rotor3<T>& r) {
        return geometric_product(m, detail::promote(r));
    }

    template<typename T>
    /**
     * @brief Geometric product of a rotor and a multivector.
     * @param r The rotor operand.
     * @param m The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator*(const Rotor3<T>& r, const Multivector3<T>& m) {
        return geometric_product(detail::promote(r), m);
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
    constexpr Trivector3<T> operator^(const Vector3<T>& v, const Bivector3<T>& b) {
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
    constexpr Trivector3<T> operator^(const Bivector3<T>& b, const Vector3<T>& v) {
        return geometric_product(detail::promote(b), detail::promote(v)).trivector;
    }

    template<typename T>
    /**
     * @brief Wedge product of a vector and a bivector.
     * @param v The vector operand.
     * @param b The bivector operand.
     * @return The resulting trivector.
     */
    constexpr Trivector3<T> wedge_product(const Vector3<T>& v, const Bivector3<T>& b) {
        return v ^ b;
    }

    template<typename T>
    /**
     * @brief Wedge product of a bivector and a vector.
     * @param b The bivector operand.
     * @param v The vector operand.
     * @return The resulting trivector.
     */
    constexpr Trivector3<T> wedge_product(const Bivector3<T>& b, const Vector3<T>& v) {
        return b ^ v;
    }
} // namespace CliffordCore
