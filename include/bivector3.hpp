#pragma once

/**
 * @file bivector3.hpp
 * @brief Grade 2: the Bivector3 type, an oriented plane element.
 *
 * Bivector3<T> is the grade 2 element: an oriented plane with a magnitude.
 * Components are stored as (xy, xz, yz) = (e1e2, e1e3, e2e3) -- note xz,
 * not the more common zx. That single choice fixes every sign in the
 * wedge, the dual and the rotation formulas.
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

#include "detail/format.hpp"

namespace CliffordCore
{

template<typename T>
class Scalar;

template<typename T>

/**
 * @brief A class representing a bivector in 3D space.
 * @tparam T The arithmetic component type.
 */
class Bivector3
{
    static_assert(std::is_arithmetic<T>::value, "Bivector3 can only be instantiated with numeric types.");

public:
    T xy;   ///< The xy-component (e1e2).
    T xz;   ///< The xz-component (e1e3).
    T yz;   ///< The yz-component (e2e3).

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
     * @brief Computes the magnitude of the bivector.
     * @return The magnitude as a scalar.
     */
    constexpr Scalar<T> magnitude() const {                                                                         
        return Scalar<T>(std::sqrt(xy * xy + xz * xz + yz * yz));
    }

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

    /**
     * @brief Adds another bivector in place.
     * @param other The bivector to add.
     * @return A reference to this bivector.
     */
    constexpr Bivector3& operator+=(const Bivector3& other) {
        xy += other.xy; xz += other.xz; yz += other.yz;
        return *this;
    }

    /**
     * @brief Subtracts another bivector in place.
     * @param other The bivector to subtract.
     * @return A reference to this bivector.
     */
    constexpr Bivector3& operator-=(const Bivector3& other) {
        xy -= other.xy; xz -= other.xz; yz -= other.yz;
        return *this;
    }

    /**
     * @brief Scales this bivector in place by a raw numeric value.
     * @param value The value to multiply by.
     * @return A reference to this bivector.
     */
    constexpr Bivector3& operator*=(T value) {
        xy *= value; xz *= value; yz *= value;
        return *this;
    }

    /**
     * @brief Divides this bivector in place by a raw numeric value.
     * @param value The value to divide by.
     * @return A reference to this bivector.
     */
    constexpr Bivector3& operator/=(T value) {
        xy /= value; xz /= value; yz /= value;
        return *this;
    }

    /**
     * @brief Scales this bivector in place by a scalar.
     * @param scalar The scalar to multiply by.
     * @return A reference to this bivector.
     */
    constexpr Bivector3& operator*=(const Scalar<T>& scalar);

    /**
     * @brief Divides this bivector in place by a scalar.
     * @param scalar The scalar to divide by.
     * @return A reference to this bivector.
     */
    constexpr Bivector3& operator/=(const Scalar<T>& scalar);

    /**
     * @brief Returns a string representation of the bivector.
     * @return A string representing the bivector.
     */
    std::string to_string() const {
        return detail::format_component(xy) + "*e12 + "
             + detail::format_component(xz) + "*e13 + "
             + detail::format_component(yz) + "*e23";
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

    template<typename T>
    constexpr Bivector3<T>& Bivector3<T>::operator*=(const Scalar<T>& scalar) {
        return *this *= scalar.value;
    }

    template<typename T>
    constexpr Bivector3<T>& Bivector3<T>::operator/=(const Scalar<T>& scalar) {
        return *this /= scalar.value;
    }

    template<typename T>
    /**
     * @brief Scales a bivector with the raw numeric value on the left.
     * @param value The value to multiply by.
     * @param b The bivector to scale.
     * @return The resulting bivector.
     */
    constexpr Bivector3<T> operator*(T value, const Bivector3<T>& b) {
        return Bivector3<T>(value * b.xy, value * b.xz, value * b.yz);
    }
} // namespace CliffordCore