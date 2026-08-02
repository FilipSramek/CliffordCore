#pragma once

#include <type_traits>

namespace CliffordCore
{

template<typename T>
class Scalar;

template<typename T>

/**
 * @brief A class representing a bivector in 3D space.
 * @property xy The xy-component of the bivector.
 * @property xz The xz-component of the bivector.
 * @property yz The yz-component of the bivector.
 */
class Bivector3
{
    static_assert(std::is_arithmetic<T>::value, "Bivector3 can only be instantiated with numeric types.");

public:
    T xy, xz, yz;

    /**
     * @brief Default constructor initializes the bivector components to zero.
     */
    constexpr Bivector3() : xy(0), xz(0), yz(0) {}                                             
    
    /**
     * @brief Constructor initializes the bivector components to the provided values.
     * @param xy_val The xy-component of the bivector.
     * @param xz_val The xz-component of the bivector.
     * @param yz_val The yz-component of the bivector.
     */
    constexpr Bivector3(T xy_val, T xz_val, T yz_val) : xy(xy_val), xz(xz_val), yz(yz_val) {}   

    /**
     * @brief Addition operator overload.
     * @param other The bivector to add.
     * @return The resulting bivector.
     */
    constexpr Bivector3 operator+(const Bivector3& other) const {                               
        return Bivector3(xy + other.xy, xz + other.xz, yz + other.yz);
    }

    /**
     * @brief Subtraction operator overload.
     * @param other The bivector to subtract.
     * @return The resulting bivector.
     */
    constexpr Bivector3 operator-(const Bivector3& other) const {                             
        return Bivector3(xy - other.xy, xz - other.xz, yz - other.yz);
    }

    /**
     * @brief Scalar multiplication operator overload.
     * @param scalar The scalar to multiply with.
     * @return The resulting bivector.
     */
    constexpr Bivector3 operator*(const Scalar<T>& scalar) const;

    /**
     * @brief Scalar division operator overload.
     * @param scalar The scalar to divide by.
     * @return The resulting bivector.
     */
    constexpr Bivector3 operator/(const Scalar<T>& scalar) const;

    /**
     * @brief Unary negation operator overload.
     * @return The resulting bivector.
     */
    constexpr Bivector3 operator-() const {                                                     
        return Bivector3(-xy, -xz, -yz);
    }
};
} // namespace CliffordCore

#include "scalar.hpp"

namespace CliffordCore
{

template<typename T>
constexpr Bivector3<T> Bivector3<T>::operator*(const Scalar<T>& scalar) const {                              
    return Bivector3<T>(xy * scalar.value, xz * scalar.value, yz * scalar.value);
}

template<typename T>
constexpr Bivector3<T> Bivector3<T>::operator/(const Scalar<T>& scalar) const {                             
    return Bivector3<T>(xy / scalar.value, xz / scalar.value, yz / scalar.value);
}
} // namespace CliffordCore