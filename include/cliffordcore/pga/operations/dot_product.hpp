#pragma once

/**
 * @file cliffordcore/pga/operations/dot_product.hpp
 * @brief Dot (inner) product of two planes.
 *
 * The symmetric inner product of two vectors, equal to the grade 0 part of
 * their geometric product. In Cl(3,0,1) it is the dot product of the two
 * plane normals: the e0 components never contribute, because e0^2 = 0, so
 * the offsets of the planes play no part. Two planes are perpendicular iff
 * their dot product is zero, whatever their positions. Across differing
 * grades the generalisation lives in contraction.hpp.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include "../scalar.hpp"
#include "../vector.hpp"

namespace CliffordCore::PGA
{

template<typename T>
/**
 * @brief Computes the dot product (inner product) of two planes.
 * @param a The first vector.
 * @param b The second vector.
 * @return The resulting scalar: the dot product of the two normals. The e0
 *         components do not contribute.
 */
constexpr Scalar<T> operator|(const Vector<T>& a, const Vector<T>& b) {
    return Scalar<T>(a.e1 * b.e1 + a.e2 * b.e2 + a.e3 * b.e3);
}

template<typename T>
/**
 * @brief Computes the dot product of two planes.
 * @param a The first vector.
 * @param b The second vector.
 * @return The resulting scalar from the dot product of vectors a and b.
 */
constexpr Scalar<T> dot_product(const Vector<T>& a, const Vector<T>& b) {
    return a | b;
}
} // namespace CliffordCore::PGA
