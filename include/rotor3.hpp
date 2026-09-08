#pragma once

#include <type_traits>
#include "scalar.hpp"
#include "bivector3.hpp"
#include "vector3.hpp"
#include "multivector3.hpp"

namespace CliffordCore
{
    template<typename T>
    /**
     * @brief A class representing a rotor in 3D space.
     * @property scalar The scalar component of the rotor.
     * @property bivector The bivector component of the rotor.
     */
    class Rotor3
    {
        static_assert(std::is_arithmetic<T>::value, "Rotor3 can only be instantiated with numeric types.");
    public:
        Scalar<T> scalar;
        Bivector3<T> bivector;

        /**
         * @brief Default constructor initializes all components to zero.
         */
        constexpr Rotor3() : scalar(), bivector() {}

        /**
         * @brief Constructor initializes the rotor with the provided components.
         * @param s The scalar component.
         * @param b The bivector component.
         */
        constexpr Rotor3(const Scalar<T>& s, const Bivector3<T>& b)
            : scalar(s), bivector(b) {}

        /**
         * @brief Constructor initializes the rotor with a 3D multivector.
         * @param m The multivector to initialize the rotor with.
         */
        constexpr Rotor3(const Multivector3<T>& m)
            : scalar(m.scalar), bivector(m.bivector) {}

        /**
         * @brief Addition operator overload.
         * @param other The rotor to add.
         * @return The resulting rotor.
         */
        constexpr Rotor3 operator+(const Rotor3& other) const {
            return Rotor3(
                scalar + other.scalar,
                bivector + other.bivector
            );
        }

        /**
         * @brief Subtraction operator overload.
         * @param other The rotor to subtract.
         * @return The resulting rotor.
         */
        constexpr Rotor3 operator-(const Rotor3& other) const {
            return Rotor3(
                scalar - other.scalar,
                bivector - other.bivector
            );
        }

        /**
         * @brief Multiplication operator overload.
         * @param other The scalar to multiply.
         * @return The resulting rotor.
         */
        constexpr Rotor3 operator*(const Scalar<T>& other) const {
            return Rotor3(
                scalar * other,
                bivector * other
            );
        }

        /**
         * @brief Division operator overload.
         * @param other The scalar to divide.
         * @return The resulting rotor.
         */
        constexpr Rotor3 operator/(const Scalar<T>& other) const {
            return Rotor3(
                scalar / other,
                bivector / other
            );
        }

        /**
         * @brief Adds another rotor in place.
         * @param other The rotor to add.
         * @return A reference to this rotor.
         */
        constexpr Rotor3& operator+=(const Rotor3& other) {
            *this = *this + other;
            return *this;
        }

        /**
         * @brief Subtracts another rotor in place.
         * @param other The rotor to subtract.
         * @return A reference to this rotor.
         */
        constexpr Rotor3& operator-=(const Rotor3& other) {
            *this = *this - other;
            return *this;
        }

        /**
         * @brief Scales this rotor in place by a scalar.
         * @param other The scalar to multiply by.
         * @return A reference to this rotor.
         */
        constexpr Rotor3& operator*=(const Scalar<T>& other) {
            *this = *this * other;
            return *this;
        }

        /**
         * @brief Divides this rotor in place by a scalar.
         * @param other The scalar to divide by.
         * @return A reference to this rotor.
         */
        constexpr Rotor3& operator/=(const Scalar<T>& other) {
            *this = *this / other;
            return *this;
        }
    };

    template<typename T>
    /**
     * @brief Scales a rotor with the raw numeric value on the left.
     * @param value The value to multiply by.
     * @param r The rotor to scale.
     * @return The resulting rotor.
     */
    constexpr Rotor3<T> operator*(T value, const Rotor3<T>& r) {
        return r * Scalar<T>(value);
    }
} // namespace CliffordCore