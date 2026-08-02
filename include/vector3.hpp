#pragma once

#include <type_traits>
#include "scalar.hpp"


namespace CliffordCore
{

template<typename T>
class Vector3
{
    static_assert(std::is_arithmetic<T>::value, "Vector3 can only be instantiated with numeric types.");

public:
    T x, y, z;

    constexpr Vector3() : x(0), y(0), z(0) {}                                       // Default constructor initializes components to zero
    constexpr Vector3(T x_val, T y_val, T z_val) : x(x_val), y(y_val), z(z_val) {}  // Constructor initializes components to the provided arguments

    constexpr Vector3 operator+(const Vector3& other) const {                       // Addition operator overload
        return Vector3(x + other.x, y + other.y, z + other.z);
    }

    constexpr Vector3 operator-(const Vector3& other) const {                       // Subtraction operator overload
        return Vector3(x - other.x, y - other.y, z - other.z);
    }

    constexpr Vector3 operator*(const Scalar<T>& scalar) const {                    // Scalar multiplication operator overload
        return Vector3(x * scalar.value, y * scalar.value, z * scalar.value);
    }

    constexpr Vector3 operator/(const Scalar<T>& scalar) const {                    // Scalar division operator overload
        return Vector3(x / scalar.value, y / scalar.value, z / scalar.value);
    }

    constexpr Vector3 operator-() const {                                           // Unary negation operator overload 
        return Vector3(-x, -y, -z);
    }

};
}