#pragma once

/**
 * @file cliffordcore/pga/operations/subtraction.hpp
 * @brief Mixed-grade subtraction.
 *
 * The exact mirror of addition.hpp: any mixed-grade difference widens to a
 * Multivector via detail::promote, and only binary pairs are defined.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include "addition.hpp"

// Mixed-grade subtraction.
//
// Same rules as addition.hpp: every mixed difference widens to a Multivector,
// only binary pairs are defined, and same-grade differences live in the type
// headers. a - b is always a + (-b), which the tests assert.

namespace CliffordCore::PGA
{
    template<typename T>
    /**
     * @brief Subtracts a scalar and a vector (plane).
     * @param a The scalar operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Scalar<T>& a, const Vector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a scalar and a bivector (line).
     * @param a The scalar operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Scalar<T>& a, const Bivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a scalar and a trivector (point).
     * @param a The scalar operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Scalar<T>& a, const Trivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a scalar and a quadvector.
     * @param a The scalar operand.
     * @param b The quadvector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Scalar<T>& a, const Quadvector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a vector (plane) and a scalar.
     * @param a The vector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Vector<T>& a, const Scalar<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a vector (plane) and a bivector (line).
     * @param a The vector operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Vector<T>& a, const Bivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a vector (plane) and a trivector (point).
     * @param a The vector operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Vector<T>& a, const Trivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a vector (plane) and a quadvector.
     * @param a The vector operand.
     * @param b The quadvector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Vector<T>& a, const Quadvector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a bivector (line) and a scalar.
     * @param a The bivector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Bivector<T>& a, const Scalar<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a bivector (line) and a vector (plane).
     * @param a The bivector operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Bivector<T>& a, const Vector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a bivector (line) and a trivector (point).
     * @param a The bivector operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Bivector<T>& a, const Trivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a bivector (line) and a quadvector.
     * @param a The bivector operand.
     * @param b The quadvector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Bivector<T>& a, const Quadvector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a trivector (point) and a scalar.
     * @param a The trivector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Trivector<T>& a, const Scalar<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a trivector (point) and a vector (plane).
     * @param a The trivector operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Trivector<T>& a, const Vector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a trivector (point) and a bivector (line).
     * @param a The trivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Trivector<T>& a, const Bivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a trivector (point) and a quadvector.
     * @param a The trivector operand.
     * @param b The quadvector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Trivector<T>& a, const Quadvector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a quadvector and a scalar.
     * @param a The quadvector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Quadvector<T>& a, const Scalar<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a quadvector and a vector (plane).
     * @param a The quadvector operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Quadvector<T>& a, const Vector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a quadvector and a bivector (line).
     * @param a The quadvector operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Quadvector<T>& a, const Bivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a quadvector and a trivector (point).
     * @param a The quadvector operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Quadvector<T>& a, const Trivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a multivector and a scalar.
     * @param a The multivector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Multivector<T>& a, const Scalar<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a multivector and a vector (plane).
     * @param a The multivector operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Multivector<T>& a, const Vector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a multivector and a bivector (line).
     * @param a The multivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Multivector<T>& a, const Bivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a multivector and a trivector (point).
     * @param a The multivector operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Multivector<T>& a, const Trivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a multivector and a quadvector.
     * @param a The multivector operand.
     * @param b The quadvector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Multivector<T>& a, const Quadvector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a scalar and a multivector.
     * @param a The scalar operand.
     * @param b The multivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Scalar<T>& a, const Multivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a vector (plane) and a multivector.
     * @param a The vector operand.
     * @param b The multivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Vector<T>& a, const Multivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a bivector (line) and a multivector.
     * @param a The bivector operand.
     * @param b The multivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Bivector<T>& a, const Multivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a trivector (point) and a multivector.
     * @param a The trivector operand.
     * @param b The multivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Trivector<T>& a, const Multivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a quadvector and a multivector.
     * @param a The quadvector operand.
     * @param b The multivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Quadvector<T>& a, const Multivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }
} // namespace CliffordCore::PGA
