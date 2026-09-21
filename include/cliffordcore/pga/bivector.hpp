#pragma once

/**
 * @file cliffordcore/pga/bivector.hpp
 * @brief Grade 2: the Bivector type, which represents a line.
 *
 * Bivector<T> is the grade 2 element of Cl(3,0,1). Its six components split
 * into a Euclidean part (e12, e13, e23), which carries the line's direction,
 * and an ideal part (e01, e02, e03), which carries its moment -- the offset
 * from the origin. A line through the origin has a zero ideal part; an ideal
 * line (a line at infinity) has a zero Euclidean part. The basis order is
 * lexicographic, (e01, e02, e03, e12, e13, e23), matching Cl3's (xy, xz, yz).
 *
 * A bivector is also the Lie algebra of the motors: exp() of a bivector is a
 * motor, and log() of a motor is a bivector. Under that reading it is a
 * screw -- an axis, an angle and a pitch -- hence the alias Twist. Line is the
 * same type under its geometric name.
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
 * @brief A class representing a grade 2 element of Cl(3,0,1): a line.
 * @tparam T The arithmetic component type.
 */
class Bivector
{
    static_assert(std::is_arithmetic<T>::value, "Bivector can only be instantiated with numeric types.");

public:
    T e01;   ///< The e0e1 component (ideal: part of the moment).
    T e02;   ///< The e0e2 component (ideal: part of the moment).
    T e03;   ///< The e0e3 component (ideal: part of the moment).
    T e12;   ///< The e1e2 component (Euclidean: the z part of the direction).
    T e13;   ///< The e1e3 component (Euclidean: minus the y part of the direction).
    T e23;   ///< The e2e3 component (Euclidean: the x part of the direction).

    /**
     * @brief Default constructor initializes the bivector components to zero.
     */
    constexpr Bivector() : e01(0), e02(0), e03(0), e12(0), e13(0), e23(0) {}

    /**
     * @brief Constructor initializes the bivector components to the provided values.
     * @param e01_val The e01 component.
     * @param e02_val The e02 component.
     * @param e03_val The e03 component.
     * @param e12_val The e12 component.
     * @param e13_val The e13 component.
     * @param e23_val The e23 component.
     */
    constexpr Bivector(T e01_val, T e02_val, T e03_val, T e12_val, T e13_val, T e23_val)
        : e01(e01_val), e02(e02_val), e03(e03_val), e12(e12_val), e13(e13_val), e23(e23_val) {}

    /**
     * @brief Computes the Euclidean magnitude of the bivector: the length of the line's direction.
     * @return The magnitude as a scalar. The ideal components do not contribute, because e0^2 = 0.
     */
    constexpr Scalar<T> magnitude() const {
        return Scalar<T>(std::sqrt(e12 * e12 + e13 * e13 + e23 * e23));
    }

    /**
     * @brief Addition operator overload.
     * @param other The bivector to add.
     * @return The resulting bivector.
     */
    constexpr Bivector operator+(const Bivector& other) const {
        return Bivector(
            e01 + other.e01, e02 + other.e02, e03 + other.e03,
            e12 + other.e12, e13 + other.e13, e23 + other.e23
        );
    }

    /**
     * @brief Subtraction operator overload.
     * @param other The bivector to subtract.
     * @return The resulting bivector.
     */
    constexpr Bivector operator-(const Bivector& other) const {
        return Bivector(
            e01 - other.e01, e02 - other.e02, e03 - other.e03,
            e12 - other.e12, e13 - other.e13, e23 - other.e23
        );
    }

    /**
     * @brief Scalar multiplication operator overload.
     * @param scalar The scalar to multiply with.
     * @return The resulting bivector.
     */
    constexpr Bivector operator*(const Scalar<T>& scalar) const;

    /**
     * @brief Scalar division operator overload.
     * @param scalar The scalar to divide by.
     * @return The resulting bivector.
     */
    constexpr Bivector operator/(const Scalar<T>& scalar) const;

    /**
     * @brief Unary negation operator overload.
     * @return The resulting bivector.
     */
    constexpr Bivector operator-() const {
        return Bivector(-e01, -e02, -e03, -e12, -e13, -e23);
    }

    /**
     * @brief Adds another bivector in place.
     * @param other The bivector to add.
     * @return A reference to this bivector.
     */
    constexpr Bivector& operator+=(const Bivector& other) {
        e01 += other.e01; e02 += other.e02; e03 += other.e03;
        e12 += other.e12; e13 += other.e13; e23 += other.e23;
        return *this;
    }

    /**
     * @brief Subtracts another bivector in place.
     * @param other The bivector to subtract.
     * @return A reference to this bivector.
     */
    constexpr Bivector& operator-=(const Bivector& other) {
        e01 -= other.e01; e02 -= other.e02; e03 -= other.e03;
        e12 -= other.e12; e13 -= other.e13; e23 -= other.e23;
        return *this;
    }

    /**
     * @brief Scales this bivector in place by a raw numeric value.
     * @param value The value to multiply by.
     * @return A reference to this bivector.
     */
    constexpr Bivector& operator*=(T value) {
        e01 *= value; e02 *= value; e03 *= value;
        e12 *= value; e13 *= value; e23 *= value;
        return *this;
    }

    /**
     * @brief Divides this bivector in place by a raw numeric value.
     * @param value The value to divide by.
     * @return A reference to this bivector.
     */
    constexpr Bivector& operator/=(T value) {
        e01 /= value; e02 /= value; e03 /= value;
        e12 /= value; e13 /= value; e23 /= value;
        return *this;
    }

    /**
     * @brief Scales this bivector in place by a scalar.
     * @param scalar The scalar to multiply by.
     * @return A reference to this bivector.
     */
    constexpr Bivector& operator*=(const Scalar<T>& scalar);

    /**
     * @brief Divides this bivector in place by a scalar.
     * @param scalar The scalar to divide by.
     * @return A reference to this bivector.
     */
    constexpr Bivector& operator/=(const Scalar<T>& scalar);

    /**
     * @brief Returns a string representation of the bivector.
     * @return A string representing the bivector.
     */
    std::string to_string() const {
        return CliffordCore::detail::format_component(e01) + "*e01 + "
             + CliffordCore::detail::format_component(e02) + "*e02 + "
             + CliffordCore::detail::format_component(e03) + "*e03 + "
             + CliffordCore::detail::format_component(e12) + "*e12 + "
             + CliffordCore::detail::format_component(e13) + "*e13 + "
             + CliffordCore::detail::format_component(e23) + "*e23";
    }
};

/**
 * @brief The grade 2 element under its geometric name: a line.
 * @tparam T The arithmetic component type.
 */
template<typename T>
using Line = Bivector<T>;

/**
 * @brief The grade 2 element read as a screw: the argument of exp() and the result of log().
 * @tparam T The arithmetic component type.
 */
template<typename T>
using Twist = Bivector<T>;

} // namespace CliffordCore::PGA

#include "scalar.hpp"

namespace CliffordCore::PGA
{

    template<typename T>
    constexpr Bivector<T> Bivector<T>::operator*(const Scalar<T>& scalar) const {
        return Bivector<T>(
            e01 * scalar.value, e02 * scalar.value, e03 * scalar.value,
            e12 * scalar.value, e13 * scalar.value, e23 * scalar.value
        );
    }

    template<typename T>
    constexpr Bivector<T> Bivector<T>::operator/(const Scalar<T>& scalar) const {
        return Bivector<T>(
            e01 / scalar.value, e02 / scalar.value, e03 / scalar.value,
            e12 / scalar.value, e13 / scalar.value, e23 / scalar.value
        );
    }

    template<typename T>
    constexpr Bivector<T>& Bivector<T>::operator*=(const Scalar<T>& scalar) {
        return *this *= scalar.value;
    }

    template<typename T>
    constexpr Bivector<T>& Bivector<T>::operator/=(const Scalar<T>& scalar) {
        return *this /= scalar.value;
    }

    template<typename T>
    /**
     * @brief Scales a bivector with the raw numeric value on the left.
     * @param value The value to multiply by.
     * @param b The bivector to scale.
     * @return The resulting bivector.
     */
    constexpr Bivector<T> operator*(T value, const Bivector<T>& b) {
        return Bivector<T>(
            value * b.e01, value * b.e02, value * b.e03,
            value * b.e12, value * b.e13, value * b.e23
        );
    }
} // namespace CliffordCore::PGA
