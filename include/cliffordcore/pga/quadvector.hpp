#pragma once

/**
 * @file cliffordcore/pga/quadvector.hpp
 * @brief Grade 4: the Quadvector pseudoscalar of Cl(3,0,1).
 *
 * Quadvector<T> is the grade 4 pseudoscalar e0123. Its sole member is
 * e0123. Because it contains e0, it squares to zero: e0123^2 = 0. That is
 * what makes Cl(3,0,1) different from Cl(3,0) in practice -- multiplying
 * by the pseudoscalar is not invertible, so the dual is a complement rather
 * than a product with e0123 (see dual.hpp), and the quadvector itself has
 * no inverse and no Euclidean norm. Two quadvectors multiply to the scalar
 * zero, always.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
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
 * @brief A class representing the grade 4 pseudoscalar of Cl(3,0,1).
 * @tparam T The arithmetic component type.
 */
class Quadvector
{
    static_assert(std::is_arithmetic<T>::value, "Quadvector can only be instantiated with numeric types.");

public:
    T e0123;   ///< The pseudoscalar component (e0e1e2e3).

    /**
     * @brief Default constructor initializes the quadvector to zero.
     */
    constexpr Quadvector() : e0123(0) {}

    /**
     * @brief Constructor initializes the quadvector with the provided value.
     * @param value The value to initialize the quadvector with.
     */
    constexpr Quadvector(const T& value) : e0123(value) {}

    /**
     * @brief Addition operator overload.
     * @param other The quadvector to add.
     * @return The resulting quadvector.
     */
    constexpr Quadvector operator+(const Quadvector& other) const {
        return Quadvector(e0123 + other.e0123);
    }

    /**
     * @brief Subtraction operator overload.
     * @param other The quadvector to subtract.
     * @return The resulting quadvector.
     */
    constexpr Quadvector operator-(const Quadvector& other) const {
        return Quadvector(e0123 - other.e0123);
    }

    /**
     * @brief Geometric product of two quadvectors.
     * @param other The quadvector to multiply with.
     * @return The scalar zero. e0123 contains e0 and e0^2 = 0, so the product
     *         of two pseudoscalars vanishes identically. Defined so that every
     *         operand pair has an operator*, as the library promises.
     */
    constexpr Scalar<T> operator*(const Quadvector& other) const;

    // Quadvector is the one PGA grade type with a single-argument converting
    // constructor, so a raw `T` could reach either operator*(Quadvector) or
    // operator*(Scalar<T>) by one user-defined conversion each -- an ambiguity.
    // These exact-match overloads take precedence and resolve it, exactly as
    // Cl3's Trivector does.

    /**
     * @brief Multiplication by a raw numeric value.
     * @param scalar The value to multiply with.
     * @return The resulting quadvector.
     */
    constexpr Quadvector operator*(T scalar) const {
        return Quadvector(e0123 * scalar);
    }

    /**
     * @brief Division by a raw numeric value.
     * @param scalar The value to divide by.
     * @return The resulting quadvector.
     */
    constexpr Quadvector operator/(T scalar) const {
        return Quadvector(e0123 / scalar);
    }

    /**
     * @brief Scalar multiplication operator overload.
     * @param scalar The scalar to multiply with.
     * @return The resulting quadvector.
     */
    constexpr Quadvector operator*(const Scalar<T>& scalar) const;

    /**
     * @brief Scalar division operator overload.
     * @param scalar The scalar to divide by.
     * @return The resulting quadvector.
     */
    constexpr Quadvector operator/(const Scalar<T>& scalar) const;

    /**
     * @brief Unary negation operator overload.
     * @return The resulting quadvector.
     */
    constexpr Quadvector operator-() const {
        return Quadvector(-e0123);
    }

    /**
     * @brief Adds another quadvector in place.
     * @param other The quadvector to add.
     * @return A reference to this quadvector.
     */
    constexpr Quadvector& operator+=(const Quadvector& other) {
        e0123 += other.e0123;
        return *this;
    }

    /**
     * @brief Subtracts another quadvector in place.
     * @param other The quadvector to subtract.
     * @return A reference to this quadvector.
     */
    constexpr Quadvector& operator-=(const Quadvector& other) {
        e0123 -= other.e0123;
        return *this;
    }

    /**
     * @brief Scales this quadvector in place by a raw numeric value.
     * @param value The value to multiply by.
     * @return A reference to this quadvector.
     */
    constexpr Quadvector& operator*=(T value) {
        e0123 *= value;
        return *this;
    }

    /**
     * @brief Divides this quadvector in place by a raw numeric value.
     * @param value The value to divide by.
     * @return A reference to this quadvector.
     */
    constexpr Quadvector& operator/=(T value) {
        e0123 /= value;
        return *this;
    }

    /**
     * @brief Scales this quadvector in place by a scalar.
     * @param scalar The scalar to multiply by.
     * @return A reference to this quadvector.
     */
    constexpr Quadvector& operator*=(const Scalar<T>& scalar);

    /**
     * @brief Divides this quadvector in place by a scalar.
     * @param scalar The scalar to divide by.
     * @return A reference to this quadvector.
     */
    constexpr Quadvector& operator/=(const Scalar<T>& scalar);

    /**
     * @brief Returns a string representation of the quadvector.
     * @return A string representing the quadvector.
     */
    std::string to_string() const {
        return CliffordCore::detail::format_component(e0123) + "*e0123";
    }
};

} // namespace CliffordCore::PGA

#include "scalar.hpp"

namespace CliffordCore::PGA
{

template<typename T>
constexpr Scalar<T> Quadvector<T>::operator*(const Quadvector<T>& other) const {
    // e0123 * e0123 contains e0 twice, and e0^2 = 0. Touch `other` so the
    // parameter is not flagged as unused.
    (void)other;
    return Scalar<T>(T(0));
}

template<typename T>
constexpr Quadvector<T> Quadvector<T>::operator*(const Scalar<T>& scalar) const {
    return Quadvector<T>(e0123 * scalar.value);
}

template<typename T>
constexpr Quadvector<T> Quadvector<T>::operator/(const Scalar<T>& scalar) const {
    return Quadvector<T>(e0123 / scalar.value);
}

template<typename T>
constexpr Quadvector<T>& Quadvector<T>::operator*=(const Scalar<T>& scalar) {
    return *this *= scalar.value;
}

template<typename T>
constexpr Quadvector<T>& Quadvector<T>::operator/=(const Scalar<T>& scalar) {
    return *this /= scalar.value;
}

template<typename T>
/**
 * @brief Scales a quadvector with the raw numeric value on the left.
 * @param value The value to multiply by.
 * @param q The quadvector to scale.
 * @return The resulting quadvector.
 */
constexpr Quadvector<T> operator*(T value, const Quadvector<T>& q) {
    return Quadvector<T>(value * q.e0123);
}
} // namespace CliffordCore::PGA
