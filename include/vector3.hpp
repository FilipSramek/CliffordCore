#pragma once

#include <type_traits>
#include "scalar.hpp"


namespace CliffordCore
{

template<typename T>

/**
 * @brief A class representing a 3D vector with components of type T.
 * @property x The x-component of the vector.
 * @property y The y-component of the vector.
 * @property z The z-component of the vector.
 */
class Vector3
{
    static_assert(std::is_arithmetic<T>::value, "Vector3 can only be instantiated with numeric types.");

public:
    T x, y, z;
    
    /**
     * @brief Default constructor initializes the vector components to zero.
     */
    constexpr Vector3() : x(0), y(0), z(0) {}
    
    /**
     * @brief Constructor initializes the vector components to the provided values.
     * @param x_val The x-component of the vector.
     * @param y_val The y-component of the vector.
     * @param z_val The z-component of the vector.
     */
    constexpr Vector3(T x_val, T y_val, T z_val) : x(x_val), y(y_val), z(z_val) {} 

    /**
     * @brief Addition operator overload.
     * @param other The vector to add.
     * @return The resulting vector.
     */
    constexpr Vector3 operator+(const Vector3& other) const {                       
        return Vector3(x + other.x, y + other.y, z + other.z);
    }

    /**
     * @brief Subtraction operator overload.
     * @param other The vector to subtract.
     * @return The resulting vector.
     */
    constexpr Vector3 operator-(const Vector3& other) const {                       
        return Vector3(x - other.x, y - other.y, z - other.z);
    }

    /**
     * @brief Scalar multiplication operator overload.
     * @param scalar The scalar to multiply with.
     * @return The resulting vector.
     */
    constexpr Vector3 operator*(const Scalar<T>& scalar) const {                    
        return Vector3(x * scalar.value, y * scalar.value, z * scalar.value);
    }

    /**
     * @brief Scalar division operator overload.
     * @param scalar The scalar to divide by.
     * @return The resulting vector.
     */
    constexpr Vector3 operator/(const Scalar<T>& scalar) const {                    
        return Vector3(x / scalar.value, y / scalar.value, z / scalar.value);
    }

    /**
     * @brief Unary negation operator overload.
     * @return The resulting vector.
     */
    constexpr Vector3 operator-() const {                                            
        return Vector3(-x, -y, -z);
    }

};
} // namespace CliffordCore