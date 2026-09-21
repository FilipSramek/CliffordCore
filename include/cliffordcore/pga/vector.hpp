#pragma once

/**
 * @file cliffordcore/pga/vector.hpp
 * @brief Grade 1: the Vector type, which represents a plane.
 *
 * Vector<T> is the grade 1 element of Cl(3,0,1): a e1 + b e2 + c e3 + d e0,
 * which represents the plane a x + b y + c z + d = 0. The e0 component is
 * the degenerate direction (e0^2 = 0) and carries the plane's offset from
 * the origin; e1, e2, e3 square to +1 and carry its normal. A plane with
 * only an e0 component is the ideal plane at infinity. The alias Plane is
 * the same type under its geometric name.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <type_traits>
#include <cmath>
#include <string>

#include "../detail/format.hpp"

namespace CliffordCore::PGA
{

template<typename T>
class Scalar;

template<typename T>

/**
 * @brief A class representing a grade 1 element of Cl(3,0,1): a plane.
 * @tparam T The arithmetic component type.
 */
class Vector
{
    static_assert(std::is_arithmetic<T>::value, "Vector can only be instantiated with numeric types.");

public:
    T e0;   ///< The degenerate component (e0^2 = 0): the plane's offset d.
    T e1;   ///< The e1 component: the normal's x part.
    T e2;   ///< The e2 component: the normal's y part.
    T e3;   ///< The e3 component: the normal's z part.

    /**
     * @brief Default constructor initializes the vector components to zero.
     */
    constexpr Vector() : e0(0), e1(0), e2(0), e3(0) {}

    /**
     * @brief Constructor initializes the vector components to the provided values.
     * @param e0_val The e0 component.
     * @param e1_val The e1 component.
     * @param e2_val The e2 component.
     * @param e3_val The e3 component.
     */
    constexpr Vector(T e0_val, T e1_val, T e2_val, T e3_val)
        : e0(e0_val), e1(e1_val), e2(e2_val), e3(e3_val) {}

    /**
     * @brief Computes the Euclidean magnitude of the vector: the length of the plane's normal.
     * @return The magnitude as a scalar. The e0 component does not contribute, because e0^2 = 0.
     */
    constexpr Scalar<T> magnitude() const {
        return Scalar<T>(std::sqrt(e1 * e1 + e2 * e2 + e3 * e3));
    }

    /**
     * @brief Addition operator overload.
     * @param other The vector to add.
     * @return The resulting vector.
     */
    constexpr Vector operator+(const Vector& other) const {
        return Vector(e0 + other.e0, e1 + other.e1, e2 + other.e2, e3 + other.e3);
    }

    /**
     * @brief Subtraction operator overload.
     * @param other The vector to subtract.
     * @return The resulting vector.
     */
    constexpr Vector operator-(const Vector& other) const {
        return Vector(e0 - other.e0, e1 - other.e1, e2 - other.e2, e3 - other.e3);
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
        return Vector(-e0, -e1, -e2, -e3);
    }

    /**
     * @brief Adds another vector in place.
     * @param other The vector to add.
     * @return A reference to this vector.
     */
    constexpr Vector& operator+=(const Vector& other) {
        e0 += other.e0; e1 += other.e1; e2 += other.e2; e3 += other.e3;
        return *this;
    }

    /**
     * @brief Subtracts another vector in place.
     * @param other The vector to subtract.
     * @return A reference to this vector.
     */
    constexpr Vector& operator-=(const Vector& other) {
        e0 -= other.e0; e1 -= other.e1; e2 -= other.e2; e3 -= other.e3;
        return *this;
    }

    /**
     * @brief Scales this vector in place by a raw numeric value.
     * @param value The value to multiply by.
     * @return A reference to this vector.
     */
    constexpr Vector& operator*=(T value) {
        e0 *= value; e1 *= value; e2 *= value; e3 *= value;
        return *this;
    }

    /**
     * @brief Divides this vector in place by a raw numeric value.
     * @param value The value to divide by.
     * @return A reference to this vector.
     */
    constexpr Vector& operator/=(T value) {
        e0 /= value; e1 /= value; e2 /= value; e3 /= value;
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
        return CliffordCore::detail::format_component(e0) + "*e0 + "
             + CliffordCore::detail::format_component(e1) + "*e1 + "
             + CliffordCore::detail::format_component(e2) + "*e2 + "
             + CliffordCore::detail::format_component(e3) + "*e3";
    }
};

/**
 * @brief The grade 1 element under its geometric name: a plane.
 * @tparam T The arithmetic component type.
 */
template<typename T>
using Plane = Vector<T>;

} // namespace CliffordCore::PGA

#include "scalar.hpp"

namespace CliffordCore::PGA
{

    template<typename T>
    constexpr Vector<T> Vector<T>::operator*(const Scalar<T>& scalar) const {
        return Vector<T>(e0 * scalar.value, e1 * scalar.value, e2 * scalar.value, e3 * scalar.value);
    }

    template<typename T>
    constexpr Vector<T> Vector<T>::operator/(const Scalar<T>& scalar) const {
        return Vector<T>(e0 / scalar.value, e1 / scalar.value, e2 / scalar.value, e3 / scalar.value);
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
        return Vector<T>(value * v.e0, value * v.e1, value * v.e2, value * v.e3);
    }
} // namespace CliffordCore::PGA
