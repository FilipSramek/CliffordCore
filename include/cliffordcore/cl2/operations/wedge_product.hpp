#pragma once

/**
 * @file cliffordcore/cl2/operations/wedge_product.hpp
 * @brief The wedge (outer) product of two vectors, as operator^ and by name.
 *
 * The oriented area two vectors span. In two dimensions that is a single
 * number times e12 -- the 2D "cross product" -- and signed_area() in
 * geometry.hpp returns it as a Scalar.
 *
 * There is no wedge beyond grade 1 with grade 1 here: a vector wedged with a
 * bivector would be grade 3, and grade 3 does not exist in Cl(2,0), so that
 * product is identically zero and deliberately left undefined.
 *
 * @author Filip Sramek
 * @version 0.3.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <type_traits>
#include "../vector.hpp"
#include "../bivector.hpp"

namespace CliffordCore::Cl2
{

template<typename T>

/**
 * @brief Computes the wedge product (exterior product) of two 2D vectors, resulting in a bivector.
 * @param a The first vector.
 * @param b The second vector.
 * @return The resulting bivector from the wedge product of vectors a and b.
 */
constexpr Bivector<T> operator^(const Vector<T>& a, const Vector<T>& b) {
    return Bivector<T>(a.x * b.y - a.y * b.x);   // xy component
}

template<typename T>
/**
 * @brief Computes the wedge product (exterior product) of two 2D vectors, resulting in a bivector.
 * @param a The first vector.
 * @param b The second vector.
 * @return The resulting bivector from the wedge product of vectors a and b.
 */
constexpr Bivector<T> wedge_product(const Vector<T>& a, const Vector<T>& b) {
    return a ^ b;
}
} // namespace CliffordCore::Cl2
