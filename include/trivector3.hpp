#pragma once

#include <type_traits>

namespace CliffordCore
{

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

};

} // namespace CliffordCore