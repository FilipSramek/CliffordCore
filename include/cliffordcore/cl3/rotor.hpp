#pragma once

/**
 * @file cliffordcore/cl3/rotor.hpp
 * @brief Grades 0 and 2: the Rotor type, which represents a rotation.
 *
 * Rotor<T> is the even subalgebra -- grades 0 and 2 -- and represents a
 * rotation. Rotors compose by multiplication without widening, and apply
 * via the sandwich product R v ~R. The constructor from Multivector is
 * explicit, because narrowing silently discards grades 1 and 3.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <type_traits>
#include <string>
#include "scalar.hpp"
#include "bivector.hpp"
#include "vector.hpp"
#include "multivector.hpp"

namespace CliffordCore::Cl3
{
    template<typename T>
    /**
     * @brief A class representing a rotor in 3D space.
     * @tparam T The arithmetic component type.
     */
    class Rotor
    {
        static_assert(std::is_arithmetic<T>::value, "Rotor can only be instantiated with numeric types.");
    public:
        Scalar<T> scalar;      ///< The grade 0 part.
        Bivector<T> bivector; ///< The grade 2 part, the rotation plane.

        /**
         * @brief Default constructor initializes all components to zero.
         */
        constexpr Rotor() : scalar(), bivector() {}

        /**
         * @brief Constructor initializes the rotor with the provided components.
         * @param s The scalar component.
         * @param b The bivector component.
         */
        constexpr Rotor(const Scalar<T>& s, const Bivector<T>& b)
            : scalar(s), bivector(b) {}

        /**
         * @brief Constructor initializes the rotor from a 3D multivector.
         * @param m The multivector to narrow.
         *
         * EXPLICIT on purpose: this discards the vector and trivector parts of
         * m. While it was implicit, `rotor + multivector` silently resolved to
         * this conversion and threw grades 1 and 3 away with no diagnostic,
         * while `multivector + rotor` kept everything. Spell the narrowing out,
         * or use to_rotor(m), which does the same thing by name.
         */
        constexpr explicit Rotor(const Multivector<T>& m)
            : scalar(m.scalar), bivector(m.bivector) {}

        /**
         * @brief Addition operator overload.
         * @param other The rotor to add.
         * @return The resulting rotor.
         */
        constexpr Rotor operator+(const Rotor& other) const {
            return Rotor(
                scalar + other.scalar,
                bivector + other.bivector
            );
        }

        /**
         * @brief Subtraction operator overload.
         * @param other The rotor to subtract.
         * @return The resulting rotor.
         */
        constexpr Rotor operator-(const Rotor& other) const {
            return Rotor(
                scalar - other.scalar,
                bivector - other.bivector
            );
        }

        /**
         * @brief Multiplication operator overload.
         * @param other The scalar to multiply.
         * @return The resulting rotor.
         */
        constexpr Rotor operator*(const Scalar<T>& other) const {
            return Rotor(
                scalar * other,
                bivector * other
            );
        }

        /**
         * @brief Division operator overload.
         * @param other The scalar to divide.
         * @return The resulting rotor.
         */
        constexpr Rotor operator/(const Scalar<T>& other) const {
            return Rotor(
                scalar / other,
                bivector / other
            );
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
            return "(" + scalar.to_string() + ") + (" + bivector.to_string() + ")";
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
} // namespace CliffordCore::Cl3