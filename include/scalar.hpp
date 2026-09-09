#pragma once

/**
 * @file scalar.hpp
 * @brief Grade 0: the Scalar type and its products with every other grade.
 *
 * Scalar<T> is the grade 0 element of Cl(3,0) -- a single number carrying
 * the algebra's identity. It multiplies with every other grade, and every
 * other grade multiplies back, so this header uses the two-phase
 * declare-then-define layout: the cross-type operators are declared here,
 * the other five headers are included after the namespace closes, and the
 * operators are defined out-of-line below.
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

namespace CliffordCore {

template<typename T> class Vector3;
template<typename T> class Bivector3;
template<typename T> class Trivector3;
template<typename T> class Multivector3;
template<typename T> class Rotor3;

template <typename T>

/**
 * @brief A class representing a scalar value of type T.
 * @tparam T The arithmetic component type.
 */
class Scalar {    
    static_assert(std::is_arithmetic<T>::value, "Scalar can only be instantiated with numeric types.");

public:
    T value;   ///< The scalar value.

    /**
     * @brief Default constructor initializes the scalar value to zero.
     */
    constexpr Scalar() : value(0) {}                            
    /**
     * @brief Constructor initializes the scalar value to the provided argument.
     * @param val The scalar value.
     */
    constexpr Scalar(T val) : value(val) {}                     

    /**
     * @brief Addition operator overload.
     * @param other The scalar to add.
     * @return The resulting scalar.
     */
    constexpr Scalar operator+(const Scalar& other) const;

    /**
     * @brief Subtraction operator overload.
     * @param other The scalar to subtract.
     * @return The resulting scalar.
     */
    constexpr Scalar operator-(const Scalar& other) const;

    /**
     * @brief Multiplication operator overload.
     * @param other The scalar to multiply with.
     * @return The resulting scalar.
     */
    constexpr Scalar operator*(const Scalar& other) const;

    /**
     * @brief Multiplication of scaler by vector operator overload.
     * @param other The vector to multiply with.
     * @return The resulting vector.
     */
    constexpr Vector3<T> operator*(const Vector3<T>& other) const;

    /**
     * @brief Multiplication of scaler by bivector operator overload.
     * @param other The bivector to multiply with.
     * @return The resulting bivector.
     */
    constexpr Bivector3<T> operator*(const Bivector3<T>& other) const;

    /**
     * @brief Multiplication of scaler by trivector operator overload.
     * @param other The trivector to multiply with.
     * @return The resulting trivector.
     */
    constexpr Trivector3<T> operator*(const Trivector3<T>& other) const;

    /**
     * @brief Multiplication of scaler by multivector operator overload.
     * @param other The multivector to multiply with.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator*(const Multivector3<T>& other) const;

    /**
     * @brief Multiplication of scaler by rotor operator overload.
     * @param other The rotor to multiply with.
     * @return The resulting rotor.
     */
    constexpr Rotor3<T> operator*(const Rotor3<T>& other) const;

    /**
     * @brief Division operator overload.
     * @param other The scalar to divide by.
     * @return The resulting scalar.
     */
    constexpr Scalar operator/(const Scalar& other) const;

    /**
     * @brief Unary negation operator overload.
     * @return The resulting scalar.
     */
    constexpr Scalar operator-() const;

    /**
     * @brief Adds another scalar in place.
     * @param other The scalar to add.
     * @return A reference to this scalar.
     */
    constexpr Scalar& operator+=(const Scalar& other) {
        value += other.value;
        return *this;
    }

    /**
     * @brief Subtracts another scalar in place.
     * @param other The scalar to subtract.
     * @return A reference to this scalar.
     */
    constexpr Scalar& operator-=(const Scalar& other) {
        value -= other.value;
        return *this;
    }

    /**
     * @brief Multiplies by another scalar in place.
     * @param other The scalar to multiply by.
     * @return A reference to this scalar.
     */
    constexpr Scalar& operator*=(const Scalar& other) {
        value *= other.value;
        return *this;
    }

    /**
     * @brief Divides by another scalar in place.
     * @param other The scalar to divide by.
     * @return A reference to this scalar.
     */
    constexpr Scalar& operator/=(const Scalar& other) {
        value /= other.value;
        return *this;
    }

    /**
     * @brief Returns a string representation of the scalar.
     * @return A string representing the scalar.
     */
    std::string to_string() const {
        return detail::format_component(value);
    }
};
} // namespace CliffordCore

#include "vector3.hpp"
#include "bivector3.hpp"
#include "trivector3.hpp"
#include "multivector3.hpp"
#include "rotor3.hpp"

namespace CliffordCore {

template<typename T>
constexpr Scalar<T> Scalar<T>::operator+(const Scalar& other) const {     
    return Scalar(value + other.value);
}

template<typename T>
constexpr Scalar<T> Scalar<T>::operator-(const Scalar& other) const {     
    return Scalar(value - other.value);
}

template<typename T>
constexpr Scalar<T> Scalar<T>::operator*(const Scalar& other) const {     
    return Scalar(value * other.value);
}

template<typename T>
constexpr Vector3<T> Scalar<T>::operator*(const Vector3<T>& other) const { 
    return Vector3<T>(value * other.x, value * other.y, value * other.z);
}

template<typename T>
constexpr Bivector3<T> Scalar<T>::operator*(const Bivector3<T>& other) const { 
    return Bivector3<T>(value * other.xy, value * other.xz, value * other.yz);
}

template<typename T>
constexpr Trivector3<T> Scalar<T>::operator*(const Trivector3<T>& other) const { 
    return Trivector3<T>(value * other.e123);
}

template<typename T>
constexpr Multivector3<T> Scalar<T>::operator*(const Multivector3<T>& other) const {
    // Delegate to the per-grade overloads above. Multiplying by the raw `value`
    // would put a T on the left of a class type, for which no operator exists.
    return Multivector3<T>(
        *this * other.scalar,
        *this * other.vector,
        *this * other.bivector,
        *this * other.trivector
    );
}

template<typename T>
constexpr Rotor3<T> Scalar<T>::operator*(const Rotor3<T>& other) const {
    // Same reasoning as the Multivector3 overload above.
    return Rotor3<T>(
        *this * other.scalar,
        *this * other.bivector
    );
}

template<typename T>
constexpr Scalar<T> Scalar<T>::operator/(const Scalar& other) const {     
    return Scalar(value / other.value);                     
}

template<typename T>
constexpr Scalar<T> Scalar<T>::operator-() const {                         
    return Scalar(-value);
}
} // namespace CliffordCore