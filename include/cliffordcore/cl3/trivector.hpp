#pragma once

/**
 * @file cliffordcore/cl3/trivector.hpp
 * @brief Grade 3: the Trivector pseudoscalar, an oriented volume element.
 *
 * Trivector<T> is the grade 3 pseudoscalar e123: an oriented volume. Its
 * sole member is e123, not value. Because e123^2 = -1, two trivectors
 * multiply to a Scalar rather than to another Trivector, and the inverse
 * carries a minus sign.
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
 * @brief A class representing a trivector in 3D space. Also known as a pseudoscalar, it represents the oriented volume element in 3D space.
 * @tparam T The arithmetic component type.
 */
class Trivector
{
    static_assert(std::is_arithmetic<T>::value, "Trivector can only be instantiated with numeric types.");

public:
    T e123;   ///< The pseudoscalar component (e1e2e3).

    /**
     * @brief Default constructor initializes the trivector to zero.
     */
    constexpr Trivector() : e123(0) {}
    
    /** 
     * @brief Constructor initializes the trivector with the provided value.
     * @param value The value to initialize the trivector with.
     */
    constexpr Trivector(const T& value) : e123(value) {}

    /**
     * @brief Computes the magnitude of the trivector.
     * @return The magnitude as a scalar.
     */
    constexpr Scalar<T> magnitude() const {
        return Scalar<T>(std::sqrt(e123 * e123));
    }

    /**
     * @brief Addition operator overload.
     * @param other The trivector to add.
     * @return The resulting trivector.
     */
    constexpr Trivector operator+(const Trivector& other) const {
        return Trivector(e123 + other.e123);
    }

    /**
     * @brief Subtraction operator overload.
     * @param other The trivector to subtract.
     * @return The resulting trivector.
     */
    constexpr Trivector operator-(const Trivector& other) const {
        return Trivector(e123 - other.e123);
    }

    /**
     * @brief Geometric product of two trivectors.
     * @param other The trivector to multiply with.
     * @return The resulting scalar. The pseudoscalar squares to -1, so
     *         (a e123)(b e123) = -ab, a pure scalar rather than a trivector.
     */
    constexpr Scalar<T> operator*(const Trivector& other) const;

    /**
     * @brief Geometric division of two trivectors.
     * @param other The trivector to divide by.
     * @return The resulting scalar. This is a * inverse(other), and the two
     *         minus signs cancel, so the result is simply a/b.
     */
    constexpr Scalar<T> operator/(const Trivector& other) const;

    // Trivector is the one grade type with a single-argument converting
    // constructor, so a raw `T` could reach either operator*(Trivector) or
    // operator*(Scalar<T>) by one user-defined conversion each -- an ambiguity.
    // These exact-match overloads take precedence and resolve it.

    /**
     * @brief Multiplication by a raw numeric value.
     * @param scalar The value to multiply with.
     * @return The resulting trivector.
     */
    constexpr Trivector operator*(T scalar) const {
        return Trivector(e123 * scalar);
    }

    /**
     * @brief Division by a raw numeric value.
     * @param scalar The value to divide by.
     * @return The resulting trivector.
     */
    constexpr Trivector operator/(T scalar) const {
        return Trivector(e123 / scalar);
    }

    /**
     * @brief Scalar multiplication operator overload.
     * @param scalar The scalar to multiply with.
     * @return The resulting trivector.
     */
    constexpr Trivector operator*(const Scalar<T>& scalar) const;

    /**
     * @brief Scalar division operator overload.
     * @param scalar The scalar to divide by.
     * @return The resulting trivector.
     */
    constexpr Trivector operator/(const Scalar<T>& scalar) const;

    /**
     * @brief Unary negation operator overload.
     * @return The resulting trivector.
     */
    constexpr Trivector operator-() const {
        return Trivector(-e123);
    }

    /**
     * @brief Adds another trivector in place.
     * @param other The trivector to add.
     * @return A reference to this trivector.
     */
    constexpr Trivector& operator+=(const Trivector& other) {
        e123 += other.e123;
        return *this;
    }

    /**
     * @brief Subtracts another trivector in place.
     * @param other The trivector to subtract.
     * @return A reference to this trivector.
     */
    constexpr Trivector& operator-=(const Trivector& other) {
        e123 -= other.e123;
        return *this;
    }

    /**
     * @brief Scales this trivector in place by a raw numeric value.
     * @param value The value to multiply by.
     * @return A reference to this trivector.
     */
    constexpr Trivector& operator*=(T value) {
        e123 *= value;
        return *this;
    }

    /**
     * @brief Divides this trivector in place by a raw numeric value.
     * @param value The value to divide by.
     * @return A reference to this trivector.
     */
    constexpr Trivector& operator/=(T value) {
        e123 /= value;
        return *this;
    }

    /**
     * @brief Scales this trivector in place by a scalar.
     * @param scalar The scalar to multiply by.
     * @return A reference to this trivector.
     */
    constexpr Trivector& operator*=(const Scalar<T>& scalar);

    /**
     * @brief Divides this trivector in place by a scalar.
     * @param scalar The scalar to divide by.
     * @return A reference to this trivector.
     */
    constexpr Trivector& operator/=(const Scalar<T>& scalar);

    /**
     * @brief Returns a string representation of the trivector.
     * @return A string representing the trivector.
     */
    std::string to_string() const {
        return CliffordCore::detail::format_component(e123) + "*e123";
    }
};

} // namespace CliffordCore::Cl3

#include "scalar.hpp"

namespace CliffordCore::Cl3
{

template<typename T>
constexpr Scalar<T> Trivector<T>::operator*(const Trivector<T>& other) const {
    // e123 * e123 = -1, so the product of two pseudoscalars is a negative scalar.
    return Scalar<T>(-(e123 * other.e123));
}

template<typename T>
constexpr Scalar<T> Trivector<T>::operator/(const Trivector<T>& other) const {
    // a/b = a * inverse(b), and inverse(b) = -b/|b|^2, so the minus signs cancel.
    return Scalar<T>(e123 / other.e123);
}

template<typename T>
constexpr Trivector<T> Trivector<T>::operator*(const Scalar<T>& scalar) const {
    return Trivector<T>(e123 * scalar.value);
}

template<typename T>
constexpr Trivector<T> Trivector<T>::operator/(const Scalar<T>& scalar) const {
    return Trivector<T>(e123 / scalar.value);
}

template<typename T>
constexpr Trivector<T>& Trivector<T>::operator*=(const Scalar<T>& scalar) {
    return *this *= scalar.value;
}

template<typename T>
constexpr Trivector<T>& Trivector<T>::operator/=(const Scalar<T>& scalar) {
    return *this /= scalar.value;
}

template<typename T>
/**
 * @brief Scales a trivector with the raw numeric value on the left.
 * @param value The value to multiply by.
 * @param t The trivector to scale.
 * @return The resulting trivector.
 */
constexpr Trivector<T> operator*(T value, const Trivector<T>& t) {
    return Trivector<T>(value * t.e123);
}
} // namespace CliffordCore::Cl3