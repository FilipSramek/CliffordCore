#pragma once

/**
 * @file cliffordcore/pga/operations/wedge_product.hpp
 * @brief Wedge (outer) product of two planes: the line where they meet.
 *
 * The antisymmetric outer product of two vectors, equal to the grade 2 part
 * of their geometric product. In Cl(3,0,1) the wedge of two planes is the
 * line in which they intersect -- the outer product IS the meet in a
 * plane-based algebra -- and it vanishes iff the planes coincide. Two
 * parallel planes meet in an ideal line (Euclidean part zero). The higher
 * wedges, which continue the meet to points and to the pseudoscalar, live
 * in mixed_products.hpp.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <type_traits>
#include "../vector.hpp"
#include "../bivector.hpp"

namespace CliffordCore::PGA
{

template<typename T>
/**
 * @brief Computes the wedge product of two planes: the line where they meet.
 * @param a The first vector.
 * @param b The second vector.
 * @return The resulting bivector, with components in lexicographic order
 *         (e01, e02, e03, e12, e13, e23); each is a_i b_j - a_j b_i.
 */
constexpr Bivector<T> operator^(const Vector<T>& a, const Vector<T>& b) {
    return Bivector<T>(
        a.e0 * b.e1 - a.e1 * b.e0,
        a.e0 * b.e2 - a.e2 * b.e0,
        a.e0 * b.e3 - a.e3 * b.e0,
        a.e1 * b.e2 - a.e2 * b.e1,
        a.e1 * b.e3 - a.e3 * b.e1,
        a.e2 * b.e3 - a.e3 * b.e2
    );
}

template<typename T>
/**
 * @brief Computes the wedge product of two planes.
 * @param a The first vector.
 * @param b The second vector.
 * @return The resulting bivector: the line in which the planes intersect.
 */
constexpr Bivector<T> wedge_product(const Vector<T>& a, const Vector<T>& b) {
    return a ^ b;
}
} // namespace CliffordCore::PGA
