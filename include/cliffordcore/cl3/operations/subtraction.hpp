#pragma once

/**
 * @file cliffordcore/cl3/operations/subtraction.hpp
 * @brief Mixed-grade subtraction.
 *
 * The mirror of addition.hpp, with the same widening rule: same-grade
 * differences keep their type, mixed-grade differences become a
 * Multivector.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include "addition.hpp"

// Mixed-grade subtraction, the counterpart to addition.hpp.
//
// Same rule as addition: subtracting two different grades cannot stay in either
// grade, so every mixed difference widens to a Multivector, and only binary
// pairs are defined -- a - b - c parses as ((a - b) - c).
//
// Same-grade differences (Vector - Vector) keep their own type and live in the
// type headers. Multivector - Multivector is a member of Multivector.
//
// These reuse detail::promote from addition.hpp.

namespace CliffordCore::Cl3
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
    constexpr Multivector<T> operator-(const Scalar<T>& a, const Vector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a scalar from a vector.
     * @param a The vector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Vector<T>& a, const Scalar<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a bivector from a scalar.
     * @param a The scalar operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Scalar<T>& a, const Bivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a scalar from a bivector.
     * @param a The bivector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Bivector<T>& a, const Scalar<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a trivector from a scalar.
     * @param a The scalar operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Scalar<T>& a, const Trivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a scalar from a trivector.
     * @param a The trivector operand.
     * @param b The scalar operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Trivector<T>& a, const Scalar<T>& b) {
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
    constexpr Multivector<T> operator-(const Vector<T>& a, const Bivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a vector from a bivector.
     * @param a The bivector operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Bivector<T>& a, const Vector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a trivector from a vector.
     * @param a The vector operand.
     * @param b The trivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Vector<T>& a, const Trivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a vector from a trivector.
     * @param a The trivector operand.
     * @param b The vector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Trivector<T>& a, const Vector<T>& b) {
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
    constexpr Multivector<T> operator-(const Bivector<T>& a, const Trivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a bivector from a trivector.
     * @param a The trivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector carrying both grades.
     */
    constexpr Multivector<T> operator-(const Trivector<T>& a, const Bivector<T>& b) {
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
    constexpr Multivector<T> operator-(const Multivector<T>& a, const Scalar<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a multivector from a scalar.
     * @param a The scalar operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator-(const Scalar<T>& a, const Multivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a vector from a multivector.
     * @param a The multivector operand.
     * @param b The vector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator-(const Multivector<T>& a, const Vector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a multivector from a vector.
     * @param a The vector operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator-(const Vector<T>& a, const Multivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a bivector from a multivector.
     * @param a The multivector operand.
     * @param b The bivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator-(const Multivector<T>& a, const Bivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a multivector from a bivector.
     * @param a The bivector operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator-(const Bivector<T>& a, const Multivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a trivector from a multivector.
     * @param a The multivector operand.
     * @param b The trivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator-(const Multivector<T>& a, const Trivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }

    template<typename T>
    /**
     * @brief Subtracts a multivector from a trivector.
     * @param a The trivector operand.
     * @param b The multivector operand.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator-(const Trivector<T>& a, const Multivector<T>& b) {
        return detail::promote(a) - detail::promote(b);
    }
} // namespace CliffordCore::Cl3
