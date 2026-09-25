#pragma once

/**
 * @file cliffordcore/cl2/multivector.hpp
 * @brief All grades at once: the Multivector type.
 *
 * Multivector<T> holds all three grades at once and is the return type of
 * any product that can leave its operands' grades. Mixed-grade addition and
 * subtraction always widen to it, so chaining stays predictable.
 *
 * @author Filip Sramek
 * @version 0.3.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <type_traits>
#include <string>
#include "scalar.hpp"
#include "vector.hpp"
#include "bivector.hpp"
#include "rotor.hpp"

namespace CliffordCore::Cl2
{
    template<typename T>
    /**
     * @brief A class representing a multivector in 2D space.
     * @tparam T The arithmetic component type.
     */
    class Multivector
    {
        static_assert(std::is_arithmetic<T>::value, "Multivector can only be instantiated with numeric types.");
    public:
        Scalar<T> scalar;        ///< The grade 0 part.
        Vector<T> vector;        ///< The grade 1 part.
        Bivector<T> bivector;    ///< The grade 2 part.

        /**
         * @brief Default constructor initializes all components to zero.
         */
        constexpr Multivector() : scalar(), vector(), bivector() {}

        /**
         * @brief Constructor initializes the multivector with the provided components.
         * @param s The scalar component.
         * @param v The vector component.
         * @param b The bivector component.
         */
        constexpr Multivector(const Scalar<T>& s, const Vector<T>& v, const Bivector<T>& b)
            : scalar(s), vector(v), bivector(b) {}

        /**
         * @brief Constructor initializes the multivector with a 2D rotor.
         * @param r The rotor to initialize the multivector with.
         */
        constexpr Multivector(const Rotor<T>& r)
            : scalar(r.scalar), vector(Vector<T>(0, 0)), bivector(r.bivector) {}

        /**
         * @brief Addition operator overload.
         * @param other The multivector to add.
         * @return The resulting multivector.
         */
        constexpr Multivector operator+(const Multivector& other) const {
            return Multivector(
                scalar + other.scalar,
                vector + other.vector,
                bivector + other.bivector
            );
        }

        /**
         * @brief Subtraction operator overload.
         * @param other The multivector to subtract.
         * @return The resulting multivector.
         */
        constexpr Multivector operator-(const Multivector& other) const {
            return Multivector(
                scalar - other.scalar,
                vector - other.vector,
                bivector - other.bivector
            );
        }

        /**
         * @brief Unary negation operator overload.
         * @return The resulting multivector.
         */
        constexpr Multivector operator-() const {
            return Multivector(
                -scalar,
                -vector,
                -bivector
            );
        }

        /**
         * @brief Scalar multiplication operator overload.
         * @param s The scalar to multiply by.
         * @return The resulting multivector.
         */
        constexpr Multivector operator*(const Scalar<T>& s) const {
            return Multivector(
                scalar * s,
                vector * s,
                bivector * s
            );
        }

        /**
         * @brief Scalar division operator overload.
         * @param s The scalar to divide by.
         * @return The resulting multivector.
         */
        constexpr Multivector operator/(const Scalar<T>& s) const {
            return Multivector(
                scalar / s,
                vector / s,
                bivector / s
            );
        }

        /**
         * @brief Adds another multivector in place.
         * @param other The multivector to add.
         * @return A reference to this multivector.
         */
        constexpr Multivector& operator+=(const Multivector& other) {
            *this = *this + other;
            return *this;
        }

        /**
         * @brief Subtracts another multivector in place.
         * @param other The multivector to subtract.
         * @return A reference to this multivector.
         */
        constexpr Multivector& operator-=(const Multivector& other) {
            *this = *this - other;
            return *this;
        }

        /**
         * @brief Scales this multivector in place by a scalar.
         * @param s The scalar to multiply by.
         * @return A reference to this multivector.
         */
        constexpr Multivector& operator*=(const Scalar<T>& s) {
            *this = *this * s;
            return *this;
        }

        /**
         * @brief Divides this multivector in place by a scalar.
         * @param s The scalar to divide by.
         * @return A reference to this multivector.
         */
        constexpr Multivector& operator/=(const Scalar<T>& s) {
            *this = *this / s;
            return *this;
        }

        /**
         * @brief Returns a string representation of the multivector.
         * @return A string representing the multivector.
         */
        std::string to_string() const {
            // Parenthesised so the joining " + " is distinguishable from the
            // "+" inside each component.
            return "(" + scalar.to_string() + ") + (" + vector.to_string()
                 + ") + (" + bivector.to_string() + ")";
        }
    };


    template<typename T>
    /**
     * @brief Scales a multivector with the raw numeric value on the left.
     * @param value The value to multiply by.
     * @param m The multivector to scale.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(T value, const Multivector<T>& m) {
        return m * Scalar<T>(value);
    }
} // namespace CliffordCore::Cl2
