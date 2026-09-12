#pragma once

/**
 * @file dual.hpp
 * @brief Duality: multiplication by the pseudoscalar.
 *
 * dual(A) = A * e123, which exchanges grade k with grade 3-k. Because
 * e123^2 = -1, applying the dual twice negates -- that sign is correct, not
 * a defect. A plane and its normal axis are dual descriptions of the same
 * thing.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include "../scalar.hpp"
#include "../vector3.hpp"
#include "../bivector3.hpp"
#include "../trivector3.hpp"

// The dual maps a grade-k element to grade 3-k by multiplying with the
// pseudoscalar I = e123. Working the basis through, with bivectors stored as
// (xy, xz, yz) = (e1e2, e1e3, e2e3):
//
//   e1 I = e23    e2 I = -e13    e3 I = e12
//   e12 I = -e3   e13 I = e2     e23 I = -e1
//   1 I = e123    e123 I = -1
//
// Because I*I = -1, applying the dual twice negates rather than returning the
// original: dual(dual(a)) == -a. That is the correct behaviour, not a defect.

namespace CliffordCore
{

template<typename T>
/**
 * @brief Computes the dual of a scalar, resulting in a trivector.
 * @param s The scalar to compute the dual of.
 * @return The resulting trivector s * e123.
 */
constexpr Trivector3<T> dual(const Scalar<T>& s) {
    return Trivector3<T>(s.value);
}

template<typename T>
/**
 * @brief Computes the dual of a 3D vector, resulting in a bivector.
 * @param v The vector to compute the dual of.
 * @return The resulting bivector v * e123, mapping x to yz, y to -xz and z to xy.
 */
constexpr Bivector3<T> dual(const Vector3<T>& v) {
    return Bivector3<T>(v.z, -v.y, v.x);
}

template<typename T>
/**
 * @brief Computes the dual of a bivector, resulting in a 3D vector.
 * @param b The bivector to compute the dual of.
 * @return The resulting vector b * e123.
 */
constexpr Vector3<T> dual(const Bivector3<T>& b) {
    return Vector3<T>(-b.yz, b.xz, -b.xy);
}

template<typename T>
/**
 * @brief Computes the dual of a trivector, resulting in a scalar.
 * @param t The trivector to compute the dual of.
 * @return The resulting scalar. The pseudoscalar squares to -1, hence the sign.
 */
constexpr Scalar<T> dual(const Trivector3<T>& t) {
    return Scalar<T>(-t.e123);
}
} // namespace CliffordCore
