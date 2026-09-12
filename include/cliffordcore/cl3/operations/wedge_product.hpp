#pragma once

/**
 * @file cliffordcore/cl3/operations/wedge_product.hpp
 * @brief Wedge (outer) product of two vectors.
 *
 * The antisymmetric outer product, which raises grade: two vectors wedge to
 * a Bivector, and a vector with a bivector to a Trivector.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <type_traits>
#include "../vector.hpp"
#include "../bivector.hpp"

namespace CliffordCore::Cl3
{

template<typename T>

/**
 * @brief Computes the wedge product (exterior product) of two 3D vectors, resulting in a bivector.
 * @param a The first vector.
 * @param b The second vector.
 * @return The resulting bivector from the wedge product of vectors a and b.
 */
constexpr Bivector<T> operator^(const Vector<T>& a, const Vector<T>& b) {
    return Bivector<T>(
        a.x * b.y - a.y * b.x,  // xy component
        a.x * b.z - a.z * b.x,  // xz component
        a.y * b.z - a.z * b.y   // yz component
    );
}

template<typename T>
/**
 * @brief Computes the wedge product (exterior product) of two 3D vectors, resulting in a bivector.
 * @param a The first vector.
 * @param b The second vector.
 * @return The resulting bivector from the wedge product of vectors a and b.
 */
constexpr Bivector<T> wedge_product(const Vector<T>& a, const Vector<T>& b) {
    return a ^ b;
}
} // namespace CliffordCore::Cl3