#pragma once

#include <type_traits>
#include <cmath>

namespace CliffordCore
{

template<typename T>
class Scalar;

template<typename T>
/**
 * @brief A class representing a trivector in 3D space. Also known as a pseudoscalar, it represents the oriented volume element in 3D space.
 * @property e123 The trivector component.
 */
class Trivector3
{
    static_assert(std::is_arithmetic<T>::value, "Trivector3 can only be instantiated with numeric types.");

public:
    T e123;

    /**
     * @brief Default constructor initializes the trivector to zero.
     */
    constexpr Trivector3() : e123(0) {}
    
    /** 
     * @brief Constructor initializes the trivector with the provided value.
     * @param value The value to initialize the trivector with.
     */
    constexpr Trivector3(const T& value) : e123(value) {}

    /**
     * @brief Computes the magnitude of the trivector.
     * @return The magnitude as a scalar.
     */
    constexpr Scalar<T> magnitude() const {
        return Scalar<T>(std::sqrt(e123 * e123));
    }

    /**
     * @brief Addition operator overload.
     * @param other The trivector to add.
     * @return The resulting trivector.
     */
    constexpr Trivector3 operator+(const Trivector3& other) const {
        return Trivector3(e123 + other.e123);
    }

    /**
     * @brief Subtraction operator overload.
     * @param other The trivector to subtract.
     * @return The resulting trivector.
     */
    constexpr Trivector3 operator-(const Trivector3& other) const {
        return Trivector3(e123 - other.e123);
    }

    /**
     * @brief Multiplication operator overload.
     * @param other The trivector to multiply with.
     * @return The resulting trivector.
     */
    constexpr Trivector3 operator*(const Trivector3& other) const {
        return Trivector3(e123 * other.e123);
    }

    /** 
     * @brief Division operator overload.
     * @param other The trivector to divide by.
     * @return The resulting trivector.
     */
    constexpr Trivector3 operator/(const Trivector3& other) const {
        return Trivector3(e123 / other.e123);
    }

    // Trivector3 is the one grade type with a single-argument converting
    // constructor, so a raw `T` could reach either operator*(Trivector3) or
    // operator*(Scalar<T>) by one user-defined conversion each -- an ambiguity.
    // These exact-match overloads take precedence and resolve it.

    /**
     * @brief Multiplication by a raw numeric value.
     * @param scalar The value to multiply with.
     * @return The resulting trivector.
     */
    constexpr Trivector3 operator*(T scalar) const {
        return Trivector3(e123 * scalar);
    }

    /**
     * @brief Division by a raw numeric value.
     * @param scalar The value to divide by.
     * @return The resulting trivector.
     */
    constexpr Trivector3 operator/(T scalar) const {
        return Trivector3(e123 / scalar);
    }

    /**
     * @brief Scalar multiplication operator overload.
     * @param scalar The scalar to multiply with.
     * @return The resulting trivector.
     */
    constexpr Trivector3 operator*(const Scalar<T>& scalar) const;

    /**
     * @brief Scalar division operator overload.
     * @param scalar The scalar to divide by.
     * @return The resulting trivector.
     */
    constexpr Trivector3 operator/(const Scalar<T>& scalar) const;

    /**
     * @brief Unary negation operator overload.
     * @return The resulting trivector.
     */
    constexpr Trivector3 operator-() const {
        return Trivector3(-e123);
    }

};

} // namespace CliffordCore

#include "scalar.hpp"

namespace CliffordCore
{

template<typename T>
constexpr Trivector3<T> Trivector3<T>::operator*(const Scalar<T>& scalar) const {
    return Trivector3<T>(e123 * scalar.value);
}

template<typename T>
constexpr Trivector3<T> Trivector3<T>::operator/(const Scalar<T>& scalar) const {
    return Trivector3<T>(e123 / scalar.value);
}
} // namespace CliffordCore