#pragma once

#include "../scalar.hpp"
#include "../vector3.hpp"
#include "../bivector3.hpp"

namespace CliffordCore
{

template<typename T>
/**
 * @brief Computes the euclidian norm (magnitude) of a 3D vector.
 * @param v The vector for which to compute the norm.
 * @return The euclidian norm (magnitude) of the vector v.
 */
constexpr Scalar<T> norm(const Vector3<T>& v) {
    return Scalar<T>(std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z));
}

template<typename T>
/**
 * @brief Computes the euclidian norm (magnitude) of a bivector in 3D space.
 * @param b The bivector for which to compute the norm.
 * @return The euclidian norm (magnitude) of the bivector b.
 */
constexpr Scalar<T> norm(const Bivector3<T>& b) {
    return Scalar<T>(std::sqrt(b.xy * b.xy + b.xz * b.xz + b.yz * b.yz));  
}
} // namespace CliffordCore