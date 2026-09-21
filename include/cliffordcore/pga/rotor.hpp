#pragma once

/**
 * @file cliffordcore/pga/rotor.hpp
 * @brief Grade 0 plus the Euclidean bivectors: the Rotor type, a rotation about an axis through the origin.
 *
 * Rotor<T> is the part of the even subalgebra that has no e0 in it: a
 * scalar and the three Euclidean bivector components e12, e13, e23. It
 * represents a rotation about a line through the origin -- the same object
 * as Cl3's Rotor, with the same sign conventions -- and composes by
 * multiplication without widening. Rotations about lines that do not pass
 * through the origin, and translations, need a Motor.
 *
 * The bivector part is stored as three raw components rather than as a
 * Bivector<T>, so "the ideal part is zero" is structural rather than a
 * precondition the constructor could silently violate. grade2(rotor)
 * returns it as a Bivector when one is needed.
 *
 * The constructors from Multivector and Motor are explicit, because
 * narrowing silently discards the grades that do not fit.
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
     * @brief A class representing a rotation about an axis through the origin.
     * @tparam T The arithmetic component type.
     */
    class Rotor
    {
        static_assert(std::is_arithmetic<T>::value, "Rotor can only be instantiated with numeric types.");
    public:
        Scalar<T> scalar;   ///< The grade 0 part.
        T e12;              ///< The e1e2 component of the rotation plane.
        T e13;              ///< The e1e3 component of the rotation plane.
        T e23;              ///< The e2e3 component of the rotation plane.

        /**
         * @brief Default constructor initializes all components to zero.
         */
        constexpr Rotor() : scalar(), e12(0), e13(0), e23(0) {}

        /**
         * @brief Constructor initializes the rotor with the provided components.
         * @param s The scalar component.
         * @param e12_val The e12 component.
         * @param e13_val The e13 component.
         * @param e23_val The e23 component.
         */
        constexpr Rotor(const Scalar<T>& s, T e12_val, T e13_val, T e23_val)
            : scalar(s), e12(e12_val), e13(e13_val), e23(e23_val) {}

        /**
         * @brief Constructor narrows a multivector to its rotor part.
         * @param m The multivector to narrow.
         *
         * EXPLICIT on purpose: this discards every grade 1, 3 and 4 component
         * and the ideal half of grade 2. Spell the narrowing out, or use
         * to_rotor(m), which does the same thing by name.
         */
        constexpr explicit Rotor(const Multivector<T>& m)
            : scalar(m.scalar), e12(m.bivector.e12), e13(m.bivector.e13), e23(m.bivector.e23) {}

        /**
         * @brief Constructor narrows a motor to its rotor part.
         * @param m The motor to narrow.
         *
         * EXPLICIT on purpose: this discards the motor's translation -- the
         * ideal bivector components and e0123. Use to_rotor(m) to say so by name.
         */
        constexpr explicit Rotor(const Motor<T>& m)
            : scalar(m.scalar), e12(m.bivector.e12), e13(m.bivector.e13), e23(m.bivector.e23) {}

        /**
         * @brief Addition operator overload.
         * @param other The rotor to add.
         * @return The resulting rotor.
         */
        constexpr Rotor operator+(const Rotor& other) const {
            return Rotor(scalar + other.scalar, e12 + other.e12, e13 + other.e13, e23 + other.e23);
        }

        /**
         * @brief Subtraction operator overload.
         * @param other The rotor to subtract.
         * @return The resulting rotor.
         */
        constexpr Rotor operator-(const Rotor& other) const {
            return Rotor(scalar - other.scalar, e12 - other.e12, e13 - other.e13, e23 - other.e23);
        }

        /**
         * @brief Unary negation operator overload.
         * @return The resulting rotor, which represents the same rotation.
         */
        constexpr Rotor operator-() const {
            return Rotor(-scalar, -e12, -e13, -e23);
        }

        /**
         * @brief Scalar multiplication operator overload.
         * @param other The scalar to multiply by.
         * @return The resulting rotor.
         */
        constexpr Rotor operator*(const Scalar<T>& other) const {
            return Rotor(scalar * other, e12 * other.value, e13 * other.value, e23 * other.value);
        }

        /**
         * @brief Scalar division operator overload.
         * @param other The scalar to divide by.
         * @return The resulting rotor.
         */
        constexpr Rotor operator/(const Scalar<T>& other) const {
            return Rotor(scalar / other, e12 / other.value, e13 / other.value, e23 / other.value);
        }

        /**
         * @brief Adds another rotor in place.
         * @param other The rotor to add.
         * @return A reference to this rotor.
         */
        constexpr Rotor& operator+=(const Rotor& other) {
            *this = *this + other;
            return *this;
        }

        /**
         * @brief Subtracts another rotor in place.
         * @param other The rotor to subtract.
         * @return A reference to this rotor.
         */
        constexpr Rotor& operator-=(const Rotor& other) {
            *this = *this - other;
            return *this;
        }

        /**
         * @brief Scales this rotor in place by a scalar.
         * @param other The scalar to multiply by.
         * @return A reference to this rotor.
         */
        constexpr Rotor& operator*=(const Scalar<T>& other) {
            *this = *this * other;
            return *this;
        }

        /**
         * @brief Divides this rotor in place by a scalar.
         * @param other The scalar to divide by.
         * @return A reference to this rotor.
         */
        constexpr Rotor& operator/=(const Scalar<T>& other) {
            *this = *this / other;
            return *this;
        }

        /**
         * @brief Returns a string representation of the rotor.
         * @return A string representing the rotor.
         */
        std::string to_string() const {
            return "(" + scalar.to_string() + ") + ("
                 + CliffordCore::detail::format_component(e12) + "*e12 + "
                 + CliffordCore::detail::format_component(e13) + "*e13 + "
                 + CliffordCore::detail::format_component(e23) + "*e23)";
        }
    };

    template<typename T>
    /**
     * @brief Scales a rotor with the raw numeric value on the left.
     * @param value The value to multiply by.
     * @param r The rotor to scale.
     * @return The resulting rotor.
     */
    constexpr Rotor<T> operator*(T value, const Rotor<T>& r) {
        return r * Scalar<T>(value);
    }
} // namespace CliffordCore::PGA
