#pragma once

/**
 * @file cliffordcore/cl2/operations/dual.hpp
 * @brief Duality: multiplication by the pseudoscalar.
 *
 * dual(A) = A * e12, which exchanges grade k with grade 2-k. Because
 * e12^2 = -1, applying the dual twice negates -- that sign is correct, not
 * a defect. Grade 1 is its own dual in two dimensions, and the dual of a
 * vector is that vector turned a quarter turn counter-clockwise.
 *
 * @author Filip Sramek
 * @version 0.3.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../multivector.hpp"

// The dual maps a grade-k element to grade 2-k by multiplying with the
// pseudoscalar I = e12. Working the basis through:
//
//   1 I = e12     e1 I = e2     e2 I = -e1     e12 I = -1
//
// Because I*I = -1, applying the dual twice negates rather than returning the
// original: dual(dual(a)) == -a on every grade. That is the correct
// behaviour, not a defect, and it is the same rule as Cl(3,0).
//
// What is new in two dimensions: grade 1 maps to grade 1. dual(v) is v turned
// by +90 degrees, (x, y) -> (-y, x), which is why perp() in geometry.hpp is
// defined as this function rather than as separate arithmetic.

namespace CliffordCore::Cl2
{

template<typename T>
/**
 * @brief Computes the dual of a scalar, resulting in a bivector.
 * @param s The scalar to compute the dual of.
 * @return The resulting bivector s * e12.
 */
constexpr Bivector<T> dual(const Scalar<T>& s) {
    return Bivector<T>(s.value);
}

template<typename T>
/**
 * @brief Computes the dual of a 2D vector, resulting in another vector.
 * @param v The vector to compute the dual of.
 * @return The resulting vector v * e12, mapping (x, y) to (-y, x): the
 *         counter-clockwise quarter turn.
 */
constexpr Vector<T> dual(const Vector<T>& v) {
    return Vector<T>(-v.y, v.x);
}

template<typename T>
/**
 * @brief Computes the dual of a bivector, resulting in a scalar.
 * @param b The bivector to compute the dual of.
 * @return The resulting scalar. The pseudoscalar squares to -1, hence the sign.
 */
constexpr Scalar<T> dual(const Bivector<T>& b) {
    return Scalar<T>(-b.xy);
}

template<typename T>
/**
 * @brief Computes the dual of a multivector, grade by grade.
 * @param m The multivector to compute the dual of.
 * @return The resulting multivector m * e12: the scalar and bivector parts
 *         swap places, and the vector part turns a quarter turn.
 */
constexpr Multivector<T> dual(const Multivector<T>& m) {
    return Multivector<T>(
        dual(m.bivector),
        dual(m.vector),
        dual(m.scalar)
    );
}
} // namespace CliffordCore::Cl2
