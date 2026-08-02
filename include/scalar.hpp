#pragma once

#include <type_traits>

namespace CliffordCore {

template<typename T> class Vector3;
template<typename T> class Bivector3;
template<typename T> class Trivector3;
template<typename T> class Multivector3;
template<typename T> class Rotor3;

template <typename T>

/**
 * @brief A class representing a scalar value of type T.
 * @property value The scalar value.
 */
class Scalar {    
    static_assert(std::is_arithmetic<T>::value, "Scalar can only be instantiated with numeric types.");

public:
    T value;

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
    return Multivector3<T>(
        value * other.scalar,
        value * other.vector,
        value * other.bivector,
        value * other.trivector
    );
}

template<typename T>
constexpr Rotor3<T> Scalar<T>::operator*(const Rotor3<T>& other) const { 
    return Rotor3<T>(
        value * other.scalar,
        value * other.bivector
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