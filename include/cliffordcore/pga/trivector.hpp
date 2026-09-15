#pragma once

/**
 * @file cliffordcore/pga/trivector.hpp
 * @brief Grade 3: the Trivector type, which represents a point.
 *
 * Trivector<T> is the grade 3 element of Cl(3,0,1). Unlike Cl3's trivector
 * it has four components, not one, and it is not the pseudoscalar: e0123 is
 * (see quadvector.hpp). The e123 component is the point's weight, and the
 * three e0-containing components carry its position: the Euclidean point
 * (x, y, z) is e123 - x e023 + y e013 - z e012 (see primitives.hpp for the
 * factory that hides those signs). A trivector with e123 = 0 is an ideal
 * point, a pure direction. The alias Point is the same type under its
 * geometric name.
 *
 * Because two points share exactly two basis indices, the product of two
 * trivectors has grades 0 and 2 -- it is a Multivector, not a Scalar as it
 * is in Cl3.
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
 * @brief A class representing a grade 3 element of Cl(3,0,1): a point.
 * @tparam T The arithmetic component type.
 */
class Trivector
{
    static_assert(std::is_arithmetic<T>::value, "Trivector can only be instantiated with numeric types.");

public:
    T e012;   ///< The e0e1e2 component (ideal: minus the z coordinate times the weight).
    T e013;   ///< The e0e1e3 component (ideal: the y coordinate times the weight).
    T e023;   ///< The e0e2e3 component (ideal: minus the x coordinate times the weight).
    T e123;   ///< The e1e2e3 component: the point's weight.

    /**
     * @brief Default constructor initializes the trivector components to zero.
     */
    constexpr Trivector() : e012(0), e013(0), e023(0), e123(0) {}

    /**
     * @brief Constructor initializes the trivector components to the provided values.
     * @param e012_val The e012 component.
     * @param e013_val The e013 component.
     * @param e023_val The e023 component.
     * @param e123_val The e123 component.
     */
    constexpr Trivector(T e012_val, T e013_val, T e023_val, T e123_val)
        : e012(e012_val), e013(e013_val), e023(e023_val), e123(e123_val) {}

    /**
     * @brief Computes the Euclidean magnitude of the trivector: the absolute weight.
     * @return The magnitude as a scalar. The ideal components do not contribute, because e0^2 = 0.
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
        return Trivector(e012 + other.e012, e013 + other.e013, e023 + other.e023, e123 + other.e123);
    }

    /**
     * @brief Subtraction operator overload.
     * @param other The trivector to subtract.
     * @return The resulting trivector.
     */
    constexpr Trivector operator-(const Trivector& other) const {
        return Trivector(e012 - other.e012, e013 - other.e013, e023 - other.e023, e123 - other.e123);
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
        return Trivector(-e012, -e013, -e023, -e123);
    }

    /**
     * @brief Adds another trivector in place.
     * @param other The trivector to add.
     * @return A reference to this trivector.
     */
    constexpr Trivector& operator+=(const Trivector& other) {
        e012 += other.e012; e013 += other.e013; e023 += other.e023; e123 += other.e123;
        return *this;
    }

    /**
     * @brief Subtracts another trivector in place.
     * @param other The trivector to subtract.
     * @return A reference to this trivector.
     */
    constexpr Trivector& operator-=(const Trivector& other) {
        e012 -= other.e012; e013 -= other.e013; e023 -= other.e023; e123 -= other.e123;
        return *this;
    }

    /**
     * @brief Scales this trivector in place by a raw numeric value.
     * @param value The value to multiply by.
     * @return A reference to this trivector.
     */
    constexpr Trivector& operator*=(T value) {
        e012 *= value; e013 *= value; e023 *= value; e123 *= value;
        return *this;
    }

    /**
     * @brief Divides this trivector in place by a raw numeric value.
     * @param value The value to divide by.
     * @return A reference to this trivector.
     */
    constexpr Trivector& operator/=(T value) {
        e012 /= value; e013 /= value; e023 /= value; e123 /= value;
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
        return CliffordCore::detail::format_component(e012) + "*e012 + "
             + CliffordCore::detail::format_component(e013) + "*e013 + "
             + CliffordCore::detail::format_component(e023) + "*e023 + "
             + CliffordCore::detail::format_component(e123) + "*e123";
    }
};

/**
 * @brief The grade 3 element under its geometric name: a point.
 * @tparam T The arithmetic component type.
 */
template<typename T>
using Point = Trivector<T>;

} // namespace CliffordCore::PGA

#include "scalar.hpp"

namespace CliffordCore::PGA
{

    template<typename T>
    constexpr Trivector<T> Trivector<T>::operator*(const Scalar<T>& scalar) const {
        return Trivector<T>(e012 * scalar.value, e013 * scalar.value, e023 * scalar.value, e123 * scalar.value);
    }

    template<typename T>
    constexpr Trivector<T> Trivector<T>::operator/(const Scalar<T>& scalar) const {
        return Trivector<T>(e012 / scalar.value, e013 / scalar.value, e023 / scalar.value, e123 / scalar.value);
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
        return Trivector<T>(value * t.e012, value * t.e013, value * t.e023, value * t.e123);
    }
} // namespace CliffordCore::PGA
