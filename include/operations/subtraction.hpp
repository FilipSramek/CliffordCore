#pragma once

/**
 * @file subtraction.hpp
 * @brief Mixed-grade subtraction.
 */

#include "addition.hpp"

// Mixed-grade subtraction, the counterpart to addition.hpp.
//
// Same rule as addition: subtracting two different grades cannot stay in either
// grade, so every mixed difference widens to a Multivector3, and only binary
// pairs are defined -- a - b - c parses as ((a - b) - c).
//
// Same-grade differences (Vector3 - Vector3) keep their own type and live in the
// type headers. Multivector3 - Multivector3 is a member of Multivector3.
//
// These reuse detail::promote from addition.hpp.

namespace CliffordCore
{
    // -----------------------------------------------------------------------
    // Scalar with the higher grades
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Subtracts a vector from a scalar.
     * @param a The scalar operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator-(const Scalar<T>& a, const Vector3<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a scalar from a vector.
     * @param a The vector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator-(const Vector3<T>& a, const Scalar<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a bivector from a scalar.
     * @param a The scalar operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator-(const Scalar<T>& a, const Bivector3<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a scalar from a bivector.
     * @param a The bivector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator-(const Bivector3<T>& a, const Scalar<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a trivector from a scalar.
     * @param a The scalar operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator-(const Scalar<T>& a, const Trivector3<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a scalar from a trivector.
     * @param a The trivector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator-(const Trivector3<T>& a, const Scalar<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    // -----------------------------------------------------------------------
    // Vector with the higher grades
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Subtracts a bivector from a vector.
     * @param a The vector operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator-(const Vector3<T>& a, const Bivector3<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a vector from a bivector.
     * @param a The bivector operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator-(const Bivector3<T>& a, const Vector3<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a trivector from a vector.
     * @param a The vector operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator-(const Vector3<T>& a, const Trivector3<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a vector from a trivector.
     * @param a The trivector operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator-(const Trivector3<T>& a, const Vector3<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    // -----------------------------------------------------------------------
    // Bivector with the trivector
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Subtracts a trivector from a bivector.
     * @param a The bivector operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator-(const Bivector3<T>& a, const Trivector3<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a bivector from a trivector.
     * @param a The trivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector3<T> operator-(const Trivector3<T>& a, const Bivector3<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    // -----------------------------------------------------------------------
    // Multivector with each single grade
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Subtracts a scalar from a multivector.
     * @param a The multivector operand.
     * @param b The scalar operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator-(const Multivector3<T>& a, const Scalar<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a multivector from a scalar.
     * @param a The scalar operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator-(const Scalar<T>& a, const Multivector3<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a vector from a multivector.
     * @param a The multivector operand.
     * @param b The vector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator-(const Multivector3<T>& a, const Vector3<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a multivector from a vector.
     * @param a The vector operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator-(const Vector3<T>& a, const Multivector3<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a bivector from a multivector.
     * @param a The multivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator-(const Multivector3<T>& a, const Bivector3<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a multivector from a bivector.
     * @param a The bivector operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator-(const Bivector3<T>& a, const Multivector3<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a trivector from a multivector.
     * @param a The multivector operand.
     * @param b The trivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator-(const Multivector3<T>& a, const Trivector3<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a multivector from a trivector.
     * @param a The trivector operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator-(const Trivector3<T>& a, const Multivector3<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }
} // namespace CliffordCore
