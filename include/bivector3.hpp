#pragma once

namespace CliffordCore
{

template<typename T>
class Bivector3
{
    static_assert(std::is_arithmetic<T>::value, "Bivector3 can only be instantiated with numeric types.");

public:
    T xy, xz, yz;

    constexpr Bivector3() : xy(0), xz(0), yz(0) {}                                              // Default constructor initializes components to zero
    constexpr Bivector3(T xy_val, T xz_val, T yz_val) : xy(xy_val), xz(xz_val), yz(yz_val) {}   // Constructor initializes components to the provided arguments

    constexpr Bivector3 operator+(const Bivector3& other) const {                               // Addition operator overload
        return Bivector3(xy + other.xy, xz + other.xz, yz + other.yz);
    }

    constexpr Bivector3 operator-(const Bivector3& other) const {                               // Subtraction operator overload
        return Bivector3(xy - other.xy, xz - other.xz, yz - other.yz);
    }

    constexpr Bivector3 operator*(const Scalar<T>& scalar) const {                              // Scalar multiplication operator overload
        return Bivector3(xy * scalar.value, xz * scalar.value, yz * scalar.value);
    }

    constexpr Bivector3 operator/(const Scalar<T>& scalar) const {                              // Scalar division operator overload
        return Bivector3(xy / scalar.value, xz / scalar.value, yz / scalar.value);
    }

    constexpr Bivector3 operator-() const {                                                     // Unary negation operator overload 
        return Bivector3(-xy, -xz, -yz);
    }
};
}