#pragma once

/**
 * @file cliffordcore/pga/translator.hpp
 * @brief Grade 0 plus the ideal bivectors: the Translator type, a translation.
 *
 * Translator<T> is the part of the even subalgebra whose bivector half is
 * purely ideal: a scalar and the three components e01, e02, e03. The
 * translator by (dx, dy, dz) is 1 - (dx e01 + dy e02 + dz e03) / 2 (see
 * translator() in motor_construction.hpp), and it applies through the same
 * sandwich as every other versor. Because every ideal bivector squares to
 * zero, translators compose by simple addition of their ideal parts:
 * (s1 + a)(s2 + b) = s1 s2 + s1 b + s2 a. A translator is normalised when
 * its scalar is 1.
 *
 * The ideal part is stored as three raw components rather than as a
 * Bivector<T>, so "the Euclidean part is zero" is structural. The
 * constructors from Multivector and Motor are explicit, because narrowing
 * silently discards the rotation.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <type_traits>
#include <string>

#include "../detail/format.hpp"
#include "scalar.hpp"
#include "bivector.hpp"
#include "multivector.hpp"
#include "motor.hpp"

namespace CliffordCore::PGA
{
    template<typename T>
    /**
     * @brief A class representing a translation.
     * @tparam T The arithmetic component type.
     */
    class Translator
    {
        static_assert(std::is_arithmetic<T>::value, "Translator can only be instantiated with numeric types.");
    public:
        Scalar<T> scalar;   ///< The grade 0 part, 1 for a normalised translator.
        T e01;              ///< The e0e1 component: minus half the x displacement.
        T e02;              ///< The e0e2 component: minus half the y displacement.
        T e03;              ///< The e0e3 component: minus half the z displacement.

        /**
         * @brief Default constructor initializes all components to zero.
         */
        constexpr Translator() : scalar(), e01(0), e02(0), e03(0) {}

        /**
         * @brief Constructor initializes the translator with the provided components.
         * @param s The scalar component.
         * @param e01_val The e01 component.
         * @param e02_val The e02 component.
         * @param e03_val The e03 component.
         */
        constexpr Translator(const Scalar<T>& s, T e01_val, T e02_val, T e03_val)
            : scalar(s), e01(e01_val), e02(e02_val), e03(e03_val) {}

        /**
         * @brief Constructor narrows a multivector to its translator part.
         * @param m The multivector to narrow.
         *
         * EXPLICIT on purpose: this discards every grade 1, 3 and 4 component
         * and the Euclidean half of grade 2. Use to_translator(m) to say so by name.
         */
        constexpr explicit Translator(const Multivector<T>& m)
            : scalar(m.scalar), e01(m.bivector.e01), e02(m.bivector.e02), e03(m.bivector.e03) {}

        /**
         * @brief Constructor narrows a motor to its translator part.
         * @param m The motor to narrow.
         *
         * EXPLICIT on purpose: this discards the motor's rotation -- the
         * Euclidean bivector components and e0123. Use to_translator(m) to say so by name.
         */
        constexpr explicit Translator(const Motor<T>& m)
            : scalar(m.scalar), e01(m.bivector.e01), e02(m.bivector.e02), e03(m.bivector.e03) {}

        /**
         * @brief Addition operator overload.
         * @param other The translator to add.
         * @return The resulting translator.
         */
        constexpr Translator operator+(const Translator& other) const {
            return Translator(scalar + other.scalar, e01 + other.e01, e02 + other.e02, e03 + other.e03);
        }

        /**
         * @brief Subtraction operator overload.
         * @param other The translator to subtract.
         * @return The resulting translator.
         */
        constexpr Translator operator-(const Translator& other) const {
            return Translator(scalar - other.scalar, e01 - other.e01, e02 - other.e02, e03 - other.e03);
        }

        /**
         * @brief Unary negation operator overload.
         * @return The resulting translator, which represents the same translation.
         */
        constexpr Translator operator-() const {
            return Translator(-scalar, -e01, -e02, -e03);
        }

        /**
         * @brief Scalar multiplication operator overload.
         * @param other The scalar to multiply by.
         * @return The resulting translator.
         */
        constexpr Translator operator*(const Scalar<T>& other) const {
            return Translator(scalar * other, e01 * other.value, e02 * other.value, e03 * other.value);
        }

        /**
         * @brief Scalar division operator overload.
         * @param other The scalar to divide by.
         * @return The resulting translator.
         */
        constexpr Translator operator/(const Scalar<T>& other) const {
            return Translator(scalar / other, e01 / other.value, e02 / other.value, e03 / other.value);
        }

        /**
         * @brief Adds another translator in place.
         * @param other The translator to add.
         * @return A reference to this translator.
         */
        constexpr Translator& operator+=(const Translator& other) {
            *this = *this + other;
            return *this;
        }

        /**
         * @brief Subtracts another translator in place.
         * @param other The translator to subtract.
         * @return A reference to this translator.
         */
        constexpr Translator& operator-=(const Translator& other) {
            *this = *this - other;
            return *this;
        }

        /**
         * @brief Scales this translator in place by a scalar.
         * @param other The scalar to multiply by.
         * @return A reference to this translator.
         */
        constexpr Translator& operator*=(const Scalar<T>& other) {
            *this = *this * other;
            return *this;
        }

        /**
         * @brief Divides this translator in place by a scalar.
         * @param other The scalar to divide by.
         * @return A reference to this translator.
         */
        constexpr Translator& operator/=(const Scalar<T>& other) {
            *this = *this / other;
            return *this;
        }

        /**
         * @brief Returns a string representation of the translator.
         * @return A string representing the translator.
         */
        std::string to_string() const {
            return "(" + scalar.to_string() + ") + ("
                 + CliffordCore::detail::format_component(e01) + "*e01 + "
                 + CliffordCore::detail::format_component(e02) + "*e02 + "
                 + CliffordCore::detail::format_component(e03) + "*e03)";
        }
    };

    template<typename T>
    /**
     * @brief Scales a translator with the raw numeric value on the left.
     * @param value The value to multiply by.
     * @param t The translator to scale.
     * @return The resulting translator.
     */
    constexpr Translator<T> operator*(T value, const Translator<T>& t) {
        return t * Scalar<T>(value);
    }
} // namespace CliffordCore::PGA
