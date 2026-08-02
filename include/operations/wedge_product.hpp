#pragma once

#include <type_traits>
#include "../vector3.hpp"
#include "../bivector3.hpp"

namespace CliffordCore
{

template<typename T>

/**
 * @brief Computes the wedge product (exterior product) of two 3D vectors, resulting in a bivector.
 * @param a The first vector.
 * @param b The second vector.
 * @return The resulting bivector from the wedge product of vectors a and b.
 */
constexpr Bivector3<T> operator^(const Vector3<T>& a, const Vector3<T>& b) {
    return Bivector3<T>(
        a.x * b.y - a.y * b.x,  // xy component
        a.x * b.z - a.z * b.x,  // xz component
        a.y * b.z - a.z * b.y   // yz component
    );
}
} // namespace CliffordCore