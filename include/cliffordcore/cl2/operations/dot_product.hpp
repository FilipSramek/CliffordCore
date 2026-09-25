#pragma once

/**
 * @file cliffordcore/cl2/operations/dot_product.hpp
 * @brief The dot product of two vectors, as operator| and by name.
 *
 * The inner product of two vectors -- the grade 0 part of their geometric
 * product. operator| is the left contraction generally (see
 * contraction.hpp); for two vectors that is exactly this dot product.
 *
 * @author Filip Sramek
 * @version 0.3.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include "../scalar.hpp"
#include "../vector.hpp"

namespace CliffordCore::Cl2
{

template<typename T>
/**
 * @brief Computes the dot product (inner product) of two 2D vectors.
 * @param a The first vector.
 * @param b The second vector.
 * @return The resulting scalar from the dot product of vectors a and b.
 */
constexpr Scalar<T> operator|(const Vector<T>& a, const Vector<T>& b) {
    return Scalar<T>(a.x * b.x + a.y * b.y);
}

template<typename T>
/**
 * @brief Computes the dot product of two 2D vectors.
 * @param a The first vector.
 * @param b The second vector.
 * @return The resulting scalar from the dot product of vectors a and b.
 */
constexpr Scalar<T> dot_product(const Vector<T>& a, const Vector<T>& b) {
    return a | b;
}
} // namespace CliffordCore::Cl2
