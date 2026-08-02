#pragma once

#include "../scalar.hpp"
#include "../vector3.hpp"

namespace CliffordCore
{

template<typename T>
constexpr Scalar<T> operator|(const Vector3<T>& a, const Vector3<T>& b) {
    return Scalar<T>(a.x * b.x + a.y * b.y + a.z * b.z);
}

template<typename T>
/**
 * @brief Computes the dot product of two 3D vectors.
 * @param a The first vector.
 * @param b The second vector.
 * @return The resulting scalar from the dot product of vectors a and b.
 */
constexpr Scalar<T> dot_product(const Vector3<T>& a, const Vector3<T>& b) {
    return a | b;
}
} // namespace CliffordCore