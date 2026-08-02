#pragma once

#include <type_traits>
#include "scalar.hpp"
#include "vector3.hpp"
#include "bivector3.hpp"
#include "trivector3.hpp"

namespace CliffordCore
{
    template<typename T>
    /**
     * @brief A class representing a multivector in 3D space.
     * @property scalar The scalar component of the multivector.
     * @property vector The vector component of the multivector.
     * @property bivector The bivector component of the multivector.
     * @property trivector The trivector component of the multivector.
     */
    class Multivector3
    {
        static_assert(std::is_arithmetic<T>::value, "Multivector3 can only be instantiated with numeric types.");
    public:
        Scalar<T> scalar;
        Vector3<T> vector;
        Bivector3<T> bivector;
        Trivector3<T> trivector;

        /**
         * @brief Default constructor initializes all components to zero.
         */
        constexpr Multivector3() : scalar(), vector(), bivector(), trivector() {}

        /**
         * @brief Constructor initializes the multivector with the provided components.
         * @param s The scalar component.
         * @param v The vector component.
         * @param b The bivector component.
         * @param t The trivector component.
         */
        constexpr Multivector3(const Scalar<T>& s, const Vector3<T>& v, const Bivector3<T>& b, const Trivector3<T>& t)
            : scalar(s), vector(v), bivector(b), trivector(t) {}

        /**
         * @brief Addition operator overload.
         * @param other The multivector to add.
         * @return The resulting multivector.
         */
        constexpr Multivector3 operator+(const Multivector3& other) const {
            return Multivector3(
                scalar + other.scalar,
                vector + other.vector,
                bivector + other.bivector,
                trivector + other.trivector
            );
        }

        /**
         * @brief Subtraction operator overload.
         * @param other The multivector to subtract.
         * @return The resulting multivector.
         */
        constexpr Multivector3 operator-(const Multivector3& other) const {
            return Multivector3(
                scalar - other.scalar,
                vector - other.vector,
                bivector - other.bivector,
                trivector - other.trivector
            );
        }

        /**
         * @brief Unary negation operator overload.
         * @return The resulting multivector.
         */
        constexpr Multivector3 operator-() const {
            return Multivector3(
                -scalar,
                -vector,
                -bivector,
                -trivector
            );
        }

        /**
         * @brief Scalar multiplication operator overload.
         * @param s The scalar to multiply by.
         * @return The resulting multivector.
         */
        constexpr Multivector3 operator*(const Scalar<T>& s) const {
            return Multivector3(
                scalar * s,
                vector * s,
                bivector * s,
                trivector * s
            );
        }

        /**
         * @brief Scalar division operator overload.
         * @param s The scalar to divide by.
         * @return The resulting multivector.
         */
        constexpr Multivector3 operator/(const Scalar<T>& s) const {
            return Multivector3(
                scalar / s,
                vector / s,
                bivector / s,
                trivector / s
            );
        }
    };
} // namespace CliffordCore