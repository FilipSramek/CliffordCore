#pragma once

/**
 * @file cliffordcore/cl2/operations/norm.hpp
 * @brief Magnitudes: norm and squared_norm for every type.
 *
 * norm is the square root of squared_norm, which is the sum of the squared
 * components. squared_norm avoids the square root where the comparison does
 * not need it. The metric is Euclidean and non-degenerate, so every
 * component counts -- unlike Cl(3,0,1), nothing is left out.
 *
 * @author Filip Sramek
 * @version 0.3.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <cmath>
#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"

namespace CliffordCore::Cl2
{

template<typename T>
/**
 * @brief Computes the euclidian norm (magnitude) of a scalar.
 * @param s The scalar for which to compute the norm.
 * @return The euclidian norm (magnitude) of the scalar s.
 */
constexpr Scalar<T> norm(const Scalar<T>& s) {
    return Scalar<T>(std::abs(s.value));
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a scalar.
 * @param s The scalar for which to compute the squared norm.
 * @return The squared magnitude of the scalar s.
 */
constexpr Scalar<T> squared_norm(const Scalar<T>& s) {
    return Scalar<T>(s.value * s.value);
}

template<typename T>
/**
 * @brief Computes the euclidian norm (magnitude) of a 2D vector.
 * @param v The vector for which to compute the norm.
 * @return The euclidian norm (magnitude) of the vector v.
 */
constexpr Scalar<T> norm(const Vector<T>& v) {
    return Scalar<T>(std::sqrt(v.x * v.x + v.y * v.y));
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a 2D vector.
 * @param v The vector for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the vector v.
 */
constexpr Scalar<T> squared_norm(const Vector<T>& v) {
    return Scalar<T>(v.x * v.x + v.y * v.y);
}

template<typename T>
/**
 * @brief Computes the euclidian norm (magnitude) of a bivector in 2D space.
 * @param b The bivector for which to compute the norm.
 * @return The euclidian norm (magnitude) of the bivector b: the unsigned area.
 */
constexpr Scalar<T> norm(const Bivector<T>& b) {
    return Scalar<T>(std::abs(b.xy));
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a bivector in 2D space.
 * @param b The bivector for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the bivector b.
 */
constexpr Scalar<T> squared_norm(const Bivector<T>& b) {
    return Scalar<T>(b.xy * b.xy);
}

template<typename T>
/**
 * @brief Computes the euclidian norm (magnitude) of a multivector in 2D space.
 * @param m The multivector for which to compute the norm.
 * @return The euclidian norm (magnitude) of the multivector m.
 */
constexpr Scalar<T> norm(const Multivector<T>& m) {
    return Scalar<T>(std::sqrt(squared_norm(m).value));
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a multivector in 2D space.
 * @param m The multivector for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the multivector m, which
 *         is the grade 0 part of m * reverse(m).
 */
constexpr Scalar<T> squared_norm(const Multivector<T>& m) {
    return Scalar<T>((m.scalar.value * m.scalar.value)
        + (m.vector.x * m.vector.x) + (m.vector.y * m.vector.y)
        + (m.bivector.xy * m.bivector.xy));
}

template<typename T>
/**
 * @brief Computes the euclidian norm (magnitude) of a rotor in 2D space.
 * @param r The rotor for which to compute the norm.
 * @return The euclidian norm (magnitude) of the rotor r -- the modulus, if the
 *         rotor is read as a complex number.
 */
constexpr Scalar<T> norm(const Rotor<T>& r) {
    return Scalar<T>(std::sqrt(squared_norm(r).value));
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a rotor in 2D space.
 * @param r The rotor for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the rotor r.
 */
constexpr Scalar<T> squared_norm(const Rotor<T>& r) {
    return Scalar<T>((r.scalar.value * r.scalar.value)
        + (r.bivector.xy * r.bivector.xy));
}

} // namespace CliffordCore::Cl2
