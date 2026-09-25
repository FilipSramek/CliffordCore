#pragma once

/**
 * @file cliffordcore/cl2/bivector.hpp
 * @brief Grade 2: the Bivector, which in Cl(2,0) is also the pseudoscalar.
 *
 * Bivector<T> is the grade 2 element e12: an oriented area. In two
 * dimensions grade 2 is the top grade, so the bivector is also the
 * pseudoscalar -- there is no Trivector in this algebra, and the alias
 * Pseudoscalar<T> names the same type. Its sole member is xy, not value.
 *
 * Because e12^2 = -1, two bivectors multiply to a Scalar rather than to
 * another Bivector, and the inverse carries a minus sign. The one fact that
 * Cl(3,0) intuition gets wrong here: e12 ANTICOMMUTES with vectors, where
 * Cl(3,0)'s e123 commutes with everything.
 *
 * @author Filip Sramek
 * @version 0.3.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <type_traits>
#include <cmath>
#include <string>

#include "../detail/format.hpp"

namespace CliffordCore::Cl2
{

template<typename T>
class Scalar;

template<typename T>
/**
 * @brief A class representing a bivector in 2D space. Also the pseudoscalar, it represents the oriented area element in the plane.
 * @tparam T The arithmetic component type.
 */
class Bivector
{
    static_assert(std::is_arithmetic<T>::value, "Bivector can only be instantiated with numeric types.");

public:
    T xy;   ///< The bivector component (e1e2), which is also the pseudoscalar.

    /**
     * @brief Default constructor initializes the bivector to zero.
     */
    constexpr Bivector() : xy(0) {}

    /**
     * @brief Constructor initializes the bivector with the provided value.
     * @param value The value to initialize the bivector with.
     */
    constexpr Bivector(const T& value) : xy(value) {}

    /**
     * @brief Computes the magnitude of the bivector.
     * @return The magnitude as a scalar.
     */
    constexpr Scalar<T> magnitude() const {
        return Scalar<T>(std::sqrt(xy * xy));
    }

    /**
     * @brief Addition operator overload.
     * @param other The bivector to add.
     * @return The resulting bivector.
     */
    constexpr Bivector operator+(const Bivector& other) const {
        return Bivector(xy + other.xy);
    }

    /**
     * @brief Subtraction operator overload.
     * @param other The bivector to subtract.
     * @return The resulting bivector.
     */
    constexpr Bivector operator-(const Bivector& other) const {
        return Bivector(xy - other.xy);
    }

    /**
     * @brief Geometric product of two bivectors.
     * @param other The bivector to multiply with.
     * @return The resulting scalar. The pseudoscalar squares to -1, so
     *         (a e12)(b e12) = -ab, a pure scalar rather than a bivector.
     */
    constexpr Scalar<T> operator*(const Bivector& other) const;

    /**
     * @brief Geometric division of two bivectors.
     * @param other The bivector to divide by.
     * @return The resulting scalar. This is a * inverse(other), and the two
     *         minus signs cancel, so the result is simply a/b.
     */
    constexpr Scalar<T> operator/(const Bivector& other) const;

    // Bivector is the one grade type with a single-argument converting
    // constructor, so a raw `T` could reach either operator*(Bivector) or
    // operator*(Scalar<T>) by one user-defined conversion each -- an ambiguity.
    // These exact-match overloads take precedence and resolve it.

    /**
     * @brief Multiplication by a raw numeric value.
     * @param scalar The value to multiply with.
     * @return The resulting bivector.
     */
    constexpr Bivector operator*(T scalar) const {
        return Bivector(xy * scalar);
    }

    /**
     * @brief Division by a raw numeric value.
     * @param scalar The value to divide by.
     * @return The resulting bivector.
     */
    constexpr Bivector operator/(T scalar) const {
        return Bivector(xy / scalar);
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
        return Bivector(-xy);
    }

    /**
     * @brief Adds another bivector in place.
     * @param other The bivector to add.
     * @return A reference to this bivector.
     */
    constexpr Bivector& operator+=(const Bivector& other) {
        xy += other.xy;
        return *this;
    }

    /**
     * @brief Subtracts another bivector in place.
     * @param other The bivector to subtract.
     * @return A reference to this bivector.
     */
    constexpr Bivector& operator-=(const Bivector& other) {
        xy -= other.xy;
        return *this;
    }

    /**
     * @brief Scales this bivector in place by a raw numeric value.
     * @param value The value to multiply by.
     * @return A reference to this bivector.
     */
    constexpr Bivector& operator*=(T value) {
        xy *= value;
        return *this;
    }

    /**
     * @brief Divides this bivector in place by a raw numeric value.
     * @param value The value to divide by.
     * @return A reference to this bivector.
     */
    constexpr Bivector& operator/=(T value) {
        xy /= value;
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
        return CliffordCore::detail::format_component(xy) + "*e12";
    }
};

template<typename T>
/**
 * @brief The pseudoscalar of Cl(2,0), which is the bivector.
 *
 * Grade 2 is the top grade in two dimensions, so the area element and the
 * pseudoscalar are the same object. This alias exists so code that means "the
 * pseudoscalar" can say so.
 *
 * @tparam T The arithmetic component type.
 */
using Pseudoscalar = Bivector<T>;

} // namespace CliffordCore::Cl2

#include "scalar.hpp"

namespace CliffordCore::Cl2
{

template<typename T>
constexpr Scalar<T> Bivector<T>::operator*(const Bivector<T>& other) const {
    // e12 * e12 = -1, so the product of two pseudoscalars is a negative scalar.
    return Scalar<T>(-(xy * other.xy));
}

template<typename T>
constexpr Scalar<T> Bivector<T>::operator/(const Bivector<T>& other) const {
    // a/b = a * inverse(b), and inverse(b) = -b/|b|^2, so the minus signs cancel.
    return Scalar<T>(xy / other.xy);
}

template<typename T>
constexpr Bivector<T> Bivector<T>::operator*(const Scalar<T>& scalar) const {
    return Bivector<T>(xy * scalar.value);
}

template<typename T>
constexpr Bivector<T> Bivector<T>::operator/(const Scalar<T>& scalar) const {
    return Bivector<T>(xy / scalar.value);
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
    return Bivector<T>(value * b.xy);
}
} // namespace CliffordCore::Cl2
