#pragma once

#include <type_traits>
#include "scalar.hpp"
#include "bivector3.hpp"
#include "vector3.hpp"

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

    };
} // namespace CliffordCore