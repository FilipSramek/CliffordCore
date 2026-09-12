#pragma once

/**
 * @file addition.hpp
 * @brief Mixed-grade addition, and the rotor-flavoured conversions.
 *
 * Same-grade sums keep their own type; any mixed-grade sum widens to a
 * Multivector3 via detail::promote. Only binary pairs are defined, since
 * left-associativity extends them to chains of any length.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include "../scalar.hpp"
#include "../vector3.hpp"
#include "../bivector3.hpp"
#include "../trivector3.hpp"
#include "../multivector3.hpp"
#include "../rotor3.hpp"

// Mixed-grade addition.
//
// Adding two different grades cannot stay in either grade, so every mixed sum
// widens to a Multivector3. Only binary pairs are defined here: C++ evaluates
// a + b + c as ((a + b) + c), so chains of any length and any order follow from
// the pairs alone.
//
// Same-grade sums (Vector3 + Vector3) keep their own type and live in the type
// headers. Multivector3 + Multivector3 is a member of Multivector3.
//
// Every operator returns Multivector3, including Scalar + Bivector3 even though
// that pair is exactly a rotor. Overloads cannot differ by return type alone,
// and a uniform rule keeps chaining predictable. Use rotor_sum() or to_rotor()
// below when rotor semantics are wanted -- the same split the library already
// makes between geometric_product() and rotor_product().

namespace CliffordCore
{
    namespace detail
    {
        template<typename T>
        /**
         * @brief Widens a scalar into a multivector, zeroing the other grades.
         * @param s The scalar to widen.
         * @return The equivalent multivector.
         */
        constexpr Multivector3<T> promote(const Scalar<T>& s) {
            return Multivector3<T>(s, Vector3<T>(), Bivector3<T>(), Trivector3<T>());
        }

        template<typename T>
        /**
         * @brief Widens a vector into a multivector, zeroing the other grades.
         * @param v The vector to widen.
         * @return The equivalent multivector.
         */
        constexpr Multivector3<T> promote(const Vector3<T>& v) {
            return Multivector3<T>(Scalar<T>(), v, Bivector3<T>(), Trivector3<T>());
        }

        template<typename T>
        /**
         * @brief Widens a bivector into a multivector, zeroing the other grades.
         * @param b The bivector to widen.
         * @return The equivalent multivector.
         */
        constexpr Multivector3<T> promote(const Bivector3<T>& b) {
            return Multivector3<T>(Scalar<T>(), Vector3<T>(), b, Trivector3<T>());
        }

        template<typename T>
        /**
         * @brief Widens a trivector into a multivector, zeroing the other grades.
         * @param t The trivector to widen.
         * @return The equivalent multivector.
         */
        constexpr Multivector3<T> promote(const Trivector3<T>& t) {
            return Multivector3<T>(Scalar<T>(), Vector3<T>(), Bivector3<T>(), t);
        }

        template<typename T>
        /**
         * @brief Pass-through so multivector operands need no special case.
         * @param m The multivector.
         * @return The same multivector.
         */
        constexpr Multivector3<T> promote(const Multivector3<T>& m) {
            return m;
        }
    } // namespace detail

    // -----------------------------------------------------------------------
    // Scalar with the higher grades
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Adds a scalar and a vector.
     * @param a The scalar operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator+(const Scalar<T>& a, const Vector3<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a vector and a scalar.
     * @param a The vector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator+(const Vector3<T>& a, const Scalar<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a scalar and a bivector.
     * @param a The scalar operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator+(const Scalar<T>& a, const Bivector3<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a bivector and a scalar.
     * @param a The bivector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator+(const Bivector3<T>& a, const Scalar<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a scalar and a trivector.
     * @param a The scalar operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator+(const Scalar<T>& a, const Trivector3<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a trivector and a scalar.
     * @param a The trivector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator+(const Trivector3<T>& a, const Scalar<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    // -----------------------------------------------------------------------
    // Vector with the higher grades
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Adds a vector and a bivector.
     * @param a The vector operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator+(const Vector3<T>& a, const Bivector3<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a bivector and a vector.
     * @param a The bivector operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator+(const Bivector3<T>& a, const Vector3<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a vector and a trivector.
     * @param a The vector operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator+(const Vector3<T>& a, const Trivector3<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a trivector and a vector.
     * @param a The trivector operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator+(const Trivector3<T>& a, const Vector3<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    // -----------------------------------------------------------------------
    // Bivector with the trivector
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Adds a bivector and a trivector.
     * @param a The bivector operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator+(const Bivector3<T>& a, const Trivector3<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a trivector and a bivector.
     * @param a The trivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator+(const Trivector3<T>& a, const Bivector3<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    // -----------------------------------------------------------------------
    // Multivector with each single grade
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Adds a multivector and a scalar.
     * @param a The multivector operand.
     * @param b The scalar operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator+(const Multivector3<T>& a, const Scalar<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a scalar and a multivector.
     * @param a The scalar operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator+(const Scalar<T>& a, const Multivector3<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a multivector and a vector.
     * @param a The multivector operand.
     * @param b The vector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator+(const Multivector3<T>& a, const Vector3<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a vector and a multivector.
     * @param a The vector operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator+(const Vector3<T>& a, const Multivector3<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a multivector and a bivector.
     * @param a The multivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator+(const Multivector3<T>& a, const Bivector3<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a bivector and a multivector.
     * @param a The bivector operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator+(const Bivector3<T>& a, const Multivector3<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a multivector and a trivector.
     * @param a The multivector operand.
     * @param b The trivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator+(const Multivector3<T>& a, const Trivector3<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a trivector and a multivector.
     * @param a The trivector operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator+(const Trivector3<T>& a, const Multivector3<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    // -----------------------------------------------------------------------
    // Rotor-flavoured spellings
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Adds a scalar and a bivector as a rotor rather than a multivector.
     * @param s The scalar part.
     * @param b The bivector part.
     * @return The resulting rotor. Same components as s + b, packed as a rotor.
     */
    constexpr Rotor3<T> rotor_sum(const Scalar<T>& s, const Bivector3<T>& b) {
        return Rotor3<T>(s, b);
    }

    template<typename T>
    /**
     * @brief Narrows a multivector to a rotor.
     * @param m The multivector to narrow.
     * @return A rotor built from the scalar and bivector parts of m. The vector
     *         and trivector parts are discarded, since a rotor cannot hold them.
     */
    constexpr Rotor3<T> to_rotor(const Multivector3<T>& m) {
        return Rotor3<T>(m.scalar, m.bivector);
    }

    template<typename T>
    /**
     * @brief Widens a rotor to a multivector.
     * @param r The rotor to widen.
     * @return A multivector with the rotor's scalar and bivector parts, and zero
     *         vector and trivector parts.
     */
    constexpr Multivector3<T> to_multivector(const Rotor3<T>& r) {
        return Multivector3<T>(r.scalar, Vector3<T>(), r.bivector, Trivector3<T>());
    }
} // namespace CliffordCore
