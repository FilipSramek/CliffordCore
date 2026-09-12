#pragma once

/**
 * @file cliffordcore/cl3/operations/dot_product.hpp
 * @brief Dot (inner) product of two vectors.
 *
 * The symmetric inner product of two vectors, equal to the grade 0 part of
 * their geometric product. Across differing grades the generalisation lives
 * in contraction.hpp.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include "../scalar.hpp"
#include "../vector.hpp"

namespace CliffordCore::Cl3
{

template<typename T>
/**
 * @brief Computes the dot product (inner product) of two 3D vectors.
 * @param a The first vector.
 * @param b The second vector.
 * @return The resulting scalar from the dot product of vectors a and b.
 */
constexpr Scalar<T> operator|(const Vector<T>& a, const Vector<T>& b) {
    return Scalar<T>(a.x * b.x + a.y * b.y + a.z * b.z);
}

template<typename T>
/**
 * @brief Computes the dot product of two 3D vectors.
 * @param a The first vector.
 * @param b The second vector.
 * @return The resulting scalar from the dot product of vectors a and b.
 */
constexpr Scalar<T> dot_product(const Vector<T>& a, const Vector<T>& b) {
    return a | b;
}
} // namespace CliffordCore::Cl3