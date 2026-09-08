#pragma once

#include <type_traits>
#include <cmath>

namespace CliffordCore
{

template<typename T>
class Scalar;

template<typename T>

/**
 * @brief A class representing a 3D vector with components of type T.
 * @property x The x-component of the vector.
 * @property y The y-component of the vector.
 * @property z The z-component of the vector.
 */
class Vector3
{
    static_assert(std::is_arithmetic<T>::value, "Vector3 can only be instantiated with numeric types.");

public:
    T x, y, z;
    
    /**
     * @brief Default constructor initializes the vector components to zero.
     */
    constexpr Vector3() : x(0), y(0), z(0) {}
    
    /**
     * @brief Constructor initializes the vector components to the provided values.
     * @param x_val The x-component of the vector.
     * @param y_val The y-component of the vector.
     * @param z_val The z-component of the vector.
     */
    constexpr Vector3(T x_val, T y_val, T z_val) : x(x_val), y(y_val), z(z_val) {} 

    constexpr Scalar<T> magnitude() const {                                                                         
        return Scalar<T>(std::sqrt(x * x + y * y + z * z));
    }

    /**
     * @brief Addition operator overload.
     * @param other The vector to add.
     * @return The resulting vector.
     */
    constexpr Vector3 operator+(const Vector3& other) const {                       
        return Vector3(x + other.x, y + other.y, z + other.z);
    }

    /**
     * @brief Subtraction operator overload.
     * @param other The vector to subtract.
     * @return The resulting vector.
     */
    constexpr Vector3 operator-(const Vector3& other) const {                       
        return Vector3(x - other.x, y - other.y, z - other.z);
    }

    /**
     * @brief Scalar multiplication operator overload.
     * @param scalar The scalar to multiply with.
     * @return The resulting vector.
     */
    constexpr Vector3 operator*(const Scalar<T>& scalar) const;

    /**
     * @brief Scalar division operator overload.
     * @param scalar The scalar to divide by.
     * @return The resulting vector.
     */
    constexpr Vector3 operator/(const Scalar<T>& scalar) const;

    /**
     * @brief Unary negation operator overload.
     * @return The resulting vector.
     */
    constexpr Vector3 operator-() const {
        return Vector3(-x, -y, -z);
    }

    /**
     * @brief Adds another vector in place.
     * @param other The vector to add.
     * @return A reference to this vector.
     */
    constexpr Vector3& operator+=(const Vector3& other) {
        x += other.x; y += other.y; z += other.z;
        return *this;
    }

    /**
     * @brief Subtracts another vector in place.
     * @param other The vector to subtract.
     * @return A reference to this vector.
     */
    constexpr Vector3& operator-=(const Vector3& other) {
        x -= other.x; y -= other.y; z -= other.z;
        return *this;
    }

    /**
     * @brief Scales this vector in place by a raw numeric value.
     * @param value The value to multiply by.
     * @return A reference to this vector.
     */
    constexpr Vector3& operator*=(T value) {
        x *= value; y *= value; z *= value;
        return *this;
    }

    /**
     * @brief Divides this vector in place by a raw numeric value.
     * @param value The value to divide by.
     * @return A reference to this vector.
     */
    constexpr Vector3& operator/=(T value) {
        x /= value; y /= value; z /= value;
        return *this;
    }

    /**
     * @brief Scales this vector in place by a scalar.
     * @param scalar The scalar to multiply by.
     * @return A reference to this vector.
     */
    constexpr Vector3& operator*=(const Scalar<T>& scalar);

    /**
     * @brief Divides this vector in place by a scalar.
     * @param scalar The scalar to divide by.
     * @return A reference to this vector.
     */
    constexpr Vector3& operator/=(const Scalar<T>& scalar);

};
} // namespace CliffordCore

#include "scalar.hpp"

namespace CliffordCore
{

    template<typename T>
    constexpr Vector3<T> Vector3<T>::operator*(const Scalar<T>& scalar) const {                    
        return Vector3<T>(x * scalar.value, y * scalar.value, z * scalar.value);
    }

    template<typename T>
    constexpr Vector3<T> Vector3<T>::operator/(const Scalar<T>& scalar) const {                    
        return Vector3<T>(x / scalar.value, y / scalar.value, z / scalar.value);
    }

    template<typename T>
    constexpr Vector3<T>& Vector3<T>::operator*=(const Scalar<T>& scalar) {
        return *this *= scalar.value;
    }

    template<typename T>
    constexpr Vector3<T>& Vector3<T>::operator/=(const Scalar<T>& scalar) {
        return *this /= scalar.value;
    }

    template<typename T>
    /**
     * @brief Scales a vector with the raw numeric value on the left.
     * @param value The value to multiply by.
     * @param v The vector to scale.
     * @return The resulting vector, so that 2 * v reads the same as v * 2.
     */
    constexpr Vector3<T> operator*(T value, const Vector3<T>& v) {
        return Vector3<T>(value * v.x, value * v.y, value * v.z);
    }
} // namespace CliffordCore