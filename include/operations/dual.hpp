#pragma once

#include "../vector3.hpp"
#include "../bivector3.hpp"

namespace CliffordCore
{

template<typename T>
/**
 * @brief Computes the dual of a 3D vector, resulting in a bivector.
 * @param v The vector to compute the dual of.
 * @return The resulting bivector from the dual of vector v.
 */
constexpr Bivector3<T> dual(const Vector3<T>& v) {
    return Bivector3<T>(v.x, v.y, v.z);
}

template<typename T>
/**
 * @brief Computes the dual of a bivector, resulting in a 3D vector.
 * @param b The bivector to compute the dual of.
 * @return The resulting vector from the dual of bivector b.
 */
constexpr Vector3<T> dual(const Bivector3<T>& b) {
    return Vector3<T>(b.xy, b.xz, b.yz);
}
} // namespace CliffordCore