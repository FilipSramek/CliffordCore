#pragma once

/**
 * @file cliffordcore/cl3/vector.hpp
 * @brief Grade 1: the Vector type.
 *
 * Vector<T> is the grade 1 element: an oriented length along e1, e2, e3,
 * with e_i^2 = +1. Its product with another vector splits into the dot part
 * (a Scalar) and the wedge part (a Bivector), which is why the geometric
 * product of two vectors returns a Multivector.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <type_traits>
#include <cmath>
#include <string>

#include "../detail/format.hpp"

namespace CliffordCore::Cl3
{

template<typename T>
class Scalar;

template<typename T>

/**
 * @brief A class representing a 3D vector with components of type T.
 * @tparam T The arithmetic component type.
 */
class Vector
{
    static_assert(std::is_arithmetic<T>::value, "Vector can only be instantiated with numeric types.");

public:
    T x;   ///< The x-component of the vector.
    T y;   ///< The y-component of the vector.
    T z;   ///< The z-component of the vector.
    
    /**
     * @brief Default constructor initializes the vector components to zero.
     */
    constexpr Vector() : x(0), y(0), z(0) {}
    
    /**
     * @brief Constructor initializes the vector components to the provided values.
     * @param x_val The x-component of the vector.
     * @param y_val The y-component of the vector.
     * @param z_val The z-component of the vector.
     */
    constexpr Vector(T x_val, T y_val, T z_val) : x(x_val), y(y_val), z(z_val) {} 

    /**
     * @brief Computes the magnitude (length) of the vector.
     * @return The magnitude as a scalar.
     */
    constexpr Scalar<T> magnitude() const {                                                                         
        return Scalar<T>(std::sqrt(x * x + y * y + z * z));
    }

    /**
     * @brief Addition operator overload.
     * @param other The vector to add.
     * @return The resulting vector.
     */
    constexpr Vector operator+(const Vector& other) const {                       
        return Vector(x + other.x, y + other.y, z + other.z);
    }

    /**
     * @brief Subtraction operator overload.
     * @param other The vector to subtract.
     * @return The resulting vector.
     */
    constexpr Vector operator-(const Vector& other) const {                       
        return Vector(x - other.x, y - other.y, z - other.z);
    }

    /**
     * @brief Scalar multiplication operator overload.
     * @param scalar The scalar to multiply with.
     * @return The resulting vector.
     */
    constexpr Vector operator*(const Scalar<T>& scalar) const;

    /**
     * @brief Scalar division operator overload.
     * @param scalar The scalar to divide by.
     * @return The resulting vector.
     */
    constexpr Vector operator/(const Scalar<T>& scalar) const;

    /**
     * @brief Unary negation operator overload.
     * @return The resulting vector.
     */
    constexpr Vector operator-() const {
        return Vector(-x, -y, -z);
    }

    /**
     * @brief Adds another vector in place.
     * @param other The vector to add.
     * @return A reference to this vector.
     */
    constexpr Vector& operator+=(const Vector& other) {
        x += other.x; y += other.y; z += other.z;
        return *this;
    }

    /**
     * @brief Subtracts another vector in place.
     * @param other The vector to subtract.
     * @return A reference to this vector.
     */
    constexpr Vector& operator-=(const Vector& other) {
        x -= other.x; y -= other.y; z -= other.z;
        return *this;
    }

    /**
     * @brief Scales this vector in place by a raw numeric value.
     * @param value The value to multiply by.
     * @return A reference to this vector.
     */
    constexpr Vector& operator*=(T value) {
        x *= value; y *= value; z *= value;
        return *this;
    }

    /**
     * @brief Divides this vector in place by a raw numeric value.
     * @param value The value to divide by.
     * @return A reference to this vector.
     */
    constexpr Vector& operator/=(T value) {
        x /= value; y /= value; z /= value;
        return *this;
    }

    /**
     * @brief Scales this vector in place by a scalar.
     * @param scalar The scalar to multiply by.
     * @return A reference to this vector.
     */
    constexpr Vector& operator*=(const Scalar<T>& scalar);

    /**
     * @brief Divides this vector in place by a scalar.
     * @param scalar The scalar to divide by.
     * @return A reference to this vector.
     */
    constexpr Vector& operator/=(const Scalar<T>& scalar);

    /**
     * @brief Returns a string representation of the vector.
     * @return A string representing the vector.
     */
    std::string to_string() const {
        return CliffordCore::detail::format_component(x) + "*e1 + "
             + CliffordCore::detail::format_component(y) + "*e2 + "
             + CliffordCore::detail::format_component(z) + "*e3";
    }
};
} // namespace CliffordCore::Cl3

#include "scalar.hpp"

namespace CliffordCore::Cl3
{

    template<typename T>
    constexpr Vector<T> Vector<T>::operator*(const Scalar<T>& scalar) const {                    
        return Vector<T>(x * scalar.value, y * scalar.value, z * scalar.value);
    }

    template<typename T>
    constexpr Vector<T> Vector<T>::operator/(const Scalar<T>& scalar) const {                    
        return Vector<T>(x / scalar.value, y / scalar.value, z / scalar.value);
    }

    template<typename T>
    constexpr Vector<T>& Vector<T>::operator*=(const Scalar<T>& scalar) {
        return *this *= scalar.value;
    }

    template<typename T>
    constexpr Vector<T>& Vector<T>::operator/=(const Scalar<T>& scalar) {
        return *this /= scalar.value;
    }

    template<typename T>
    /**
     * @brief Scales a vector with the raw numeric value on the left.
     * @param value The value to multiply by.
     * @param v The vector to scale.
     * @return The resulting vector, so that 2 * v reads the same as v * 2.
     */
    constexpr Vector<T> operator*(T value, const Vector<T>& v) {
        return Vector<T>(value * v.x, value * v.y, value * v.z);
    }
} // namespace CliffordCore::Cl3