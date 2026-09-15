#pragma once

/**
 * @file cliffordcore/pga/operations/addition.hpp
 * @brief Mixed-grade addition, and the named conversions between the even types.
 *
 * Same-grade sums keep their own type; any mixed-grade sum widens to a
 * Multivector via detail::promote. Only binary pairs are defined, since
 * left-associativity extends them to chains of any length.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../quadvector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"
#include "../translator.hpp"
#include "../motor.hpp"

// Mixed-grade addition.
//
// Adding two different grades cannot stay in either grade, so every mixed sum
// widens to a Multivector. Only binary pairs are defined here: C++ evaluates
// a + b + c as ((a + b) + c), so chains of any length and any order follow from
// the pairs alone. The 20 ordered pairs of distinct single grades and the 10
// Multivector-with-a-grade pairs are all here.
//
// Same-grade sums keep their own type and live in the type headers, as do
// Rotor + Rotor, Translator + Translator and Motor + Motor. A motor plus a
// rotor or translator compiles through the implicit widening constructor
// (motor + rotor works, rotor + motor does not -- the left operand of a member
// operator never converts), and so does multivector + any even type. There is
// deliberately no Rotor + Translator: spell it to_motor(r) + to_motor(t).
//
// Every operator returns Multivector, including Scalar + Bivector even though
// that pair is exactly a motor. Overloads cannot differ by return type alone,
// and a uniform rule keeps chaining predictable. Use the to_*() conversions
// below when the narrower type is wanted.

namespace CliffordCore::PGA
{
    namespace detail
    {
        template<typename T>
        /**
         * @brief Widens a scalar into a multivector, zeroing the other grades.
         * @param s The scalar to widen.
         * @return The equivalent multivector.
         */
        constexpr Multivector<T> promote(const Scalar<T>& s) {
            return Multivector<T>(s, Vector<T>(), Bivector<T>(), Trivector<T>(), Quadvector<T>());
        }

        template<typename T>
        /**
         * @brief Widens a vector into a multivector, zeroing the other grades.
         * @param v The vector to widen.
         * @return The equivalent multivector.
         */
        constexpr Multivector<T> promote(const Vector<T>& v) {
            return Multivector<T>(Scalar<T>(), v, Bivector<T>(), Trivector<T>(), Quadvector<T>());
        }

        template<typename T>
        /**
         * @brief Widens a bivector into a multivector, zeroing the other grades.
         * @param b The bivector to widen.
         * @return The equivalent multivector.
         */
        constexpr Multivector<T> promote(const Bivector<T>& b) {
            return Multivector<T>(Scalar<T>(), Vector<T>(), b, Trivector<T>(), Quadvector<T>());
        }

        template<typename T>
        /**
         * @brief Widens a trivector into a multivector, zeroing the other grades.
         * @param t The trivector to widen.
         * @return The equivalent multivector.
         */
        constexpr Multivector<T> promote(const Trivector<T>& t) {
            return Multivector<T>(Scalar<T>(), Vector<T>(), Bivector<T>(), t, Quadvector<T>());
        }

        template<typename T>
        /**
         * @brief Widens a quadvector into a multivector, zeroing the other grades.
         * @param q The quadvector to widen.
         * @return The equivalent multivector.
         */
        constexpr Multivector<T> promote(const Quadvector<T>& q) {
            return Multivector<T>(Scalar<T>(), Vector<T>(), Bivector<T>(), Trivector<T>(), q);
        }

        template<typename T>
        /**
         * @brief Pass-through so multivector operands need no special case.
         * @param m The multivector.
         * @return The same multivector.
         */
        constexpr Multivector<T> promote(const Multivector<T>& m) {
            return m;
        }
    } // namespace detail

    // -----------------------------------------------------------------------
    // Every ordered pair of distinct grades, then Multivector with each grade
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Adds a scalar and a vector (plane).
     * @param a The scalar operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Scalar<T>& a, const Vector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a scalar and a bivector (line).
     * @param a The scalar operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Scalar<T>& a, const Bivector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a scalar and a trivector (point).
     * @param a The scalar operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Scalar<T>& a, const Trivector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a scalar and a quadvector.
     * @param a The scalar operand.
     * @param b The quadvector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Scalar<T>& a, const Quadvector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a vector (plane) and a scalar.
     * @param a The vector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Vector<T>& a, const Scalar<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a vector (plane) and a bivector (line).
     * @param a The vector operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Vector<T>& a, const Bivector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a vector (plane) and a trivector (point).
     * @param a The vector operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Vector<T>& a, const Trivector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a vector (plane) and a quadvector.
     * @param a The vector operand.
     * @param b The quadvector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Vector<T>& a, const Quadvector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a bivector (line) and a scalar.
     * @param a The bivector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Bivector<T>& a, const Scalar<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a bivector (line) and a vector (plane).
     * @param a The bivector operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Bivector<T>& a, const Vector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a bivector (line) and a trivector (point).
     * @param a The bivector operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Bivector<T>& a, const Trivector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a bivector (line) and a quadvector.
     * @param a The bivector operand.
     * @param b The quadvector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Bivector<T>& a, const Quadvector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a trivector (point) and a scalar.
     * @param a The trivector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Trivector<T>& a, const Scalar<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a trivector (point) and a vector (plane).
     * @param a The trivector operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Trivector<T>& a, const Vector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a trivector (point) and a bivector (line).
     * @param a The trivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Trivector<T>& a, const Bivector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a trivector (point) and a quadvector.
     * @param a The trivector operand.
     * @param b The quadvector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Trivector<T>& a, const Quadvector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a quadvector and a scalar.
     * @param a The quadvector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Quadvector<T>& a, const Scalar<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a quadvector and a vector (plane).
     * @param a The quadvector operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Quadvector<T>& a, const Vector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a quadvector and a bivector (line).
     * @param a The quadvector operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Quadvector<T>& a, const Bivector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a quadvector and a trivector (point).
     * @param a The quadvector operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Quadvector<T>& a, const Trivector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a multivector and a scalar.
     * @param a The multivector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Multivector<T>& a, const Scalar<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a multivector and a vector (plane).
     * @param a The multivector operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Multivector<T>& a, const Vector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a multivector and a bivector (line).
     * @param a The multivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Multivector<T>& a, const Bivector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a multivector and a trivector (point).
     * @param a The multivector operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Multivector<T>& a, const Trivector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a multivector and a quadvector.
     * @param a The multivector operand.
     * @param b The quadvector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Multivector<T>& a, const Quadvector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a scalar and a multivector.
     * @param a The scalar operand.
     * @param b The multivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Scalar<T>& a, const Multivector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a vector (plane) and a multivector.
     * @param a The vector operand.
     * @param b The multivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Vector<T>& a, const Multivector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a bivector (line) and a multivector.
     * @param a The bivector operand.
     * @param b The multivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Bivector<T>& a, const Multivector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a trivector (point) and a multivector.
     * @param a The trivector operand.
     * @param b The multivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Trivector<T>& a, const Multivector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Adds a quadvector and a multivector.
     * @param a The quadvector operand.
     * @param b The multivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator+(const Quadvector<T>& a, const Multivector<T>& b) {
        return detail::promote(a) + detail::promote(b);
    }

    // -----------------------------------------------------------------------
    // Named conversions between the even types
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Widens a rotor to a motor.
     * @param r The rotor to widen.
     * @return The equivalent motor; nothing is lost.
     */
    constexpr Motor<T> to_motor(const Rotor<T>& r) {
        return Motor<T>(r);
    }

    template<typename T>
    /**
     * @brief Widens a translator to a motor.
     * @param t The translator to widen.
     * @return The equivalent motor; nothing is lost.
     */
    constexpr Motor<T> to_motor(const Translator<T>& t) {
        return Motor<T>(t);
    }

    template<typename T>
    /**
     * @brief Narrows a multivector to a motor.
     * @param m The multivector to narrow.
     * @return A motor built from the even grades of m. Grades 1 and 3 are
     *         discarded, since a motor cannot hold them.
     */
    constexpr Motor<T> to_motor(const Multivector<T>& m) {
        return Motor<T>(m);
    }

    template<typename T>
    /**
     * @brief Narrows a multivector to a rotor.
     * @param m The multivector to narrow.
     * @return A rotor built from the scalar and the Euclidean bivector
     *         components of m. Everything else is discarded.
     */
    constexpr Rotor<T> to_rotor(const Multivector<T>& m) {
        return Rotor<T>(m);
    }

    template<typename T>
    /**
     * @brief Narrows a motor to a rotor.
     * @param m The motor to narrow.
     * @return A rotor built from the scalar and the Euclidean bivector
     *         components of m. The translation is discarded.
     */
    constexpr Rotor<T> to_rotor(const Motor<T>& m) {
        return Rotor<T>(m);
    }

    template<typename T>
    /**
     * @brief Narrows a multivector to a translator.
     * @param m The multivector to narrow.
     * @return A translator built from the scalar and the ideal bivector
     *         components of m. Everything else is discarded.
     */
    constexpr Translator<T> to_translator(const Multivector<T>& m) {
        return Translator<T>(m);
    }

    template<typename T>
    /**
     * @brief Narrows a motor to a translator.
     * @param m The motor to narrow.
     * @return A translator built from the scalar and the ideal bivector
     *         components of m. The rotation is discarded.
     */
    constexpr Translator<T> to_translator(const Motor<T>& m) {
        return Translator<T>(m);
    }

    template<typename T>
    /**
     * @brief Widens a rotor to a multivector.
     * @param r The rotor to widen.
     * @return The equivalent multivector.
     */
    constexpr Multivector<T> to_multivector(const Rotor<T>& r) {
        return Multivector<T>(r);
    }

    template<typename T>
    /**
     * @brief Widens a translator to a multivector.
     * @param t The translator to widen.
     * @return The equivalent multivector.
     */
    constexpr Multivector<T> to_multivector(const Translator<T>& t) {
        return Multivector<T>(t);
    }

    template<typename T>
    /**
     * @brief Widens a motor to a multivector.
     * @param m The motor to widen.
     * @return The equivalent multivector.
     */
    constexpr Multivector<T> to_multivector(const Motor<T>& m) {
        return Multivector<T>(m);
    }
} // namespace CliffordCore::PGA
