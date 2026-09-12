#pragma once

/**
 * @file cliffordcore/cl3/operations/norm.hpp
 * @brief Magnitudes: norm and squared_norm for every type.
 *
 * norm is the square root of squared_norm, which is the sum of the squared
 * components. squared_norm avoids the square root where the comparison does
 * not need it.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <cmath>
#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"

namespace CliffordCore::Cl3
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
 * @brief Computes the euclidian norm (magnitude) of a 3D vector.
 * @param v The vector for which to compute the norm.
 * @return The euclidian norm (magnitude) of the vector v.
 */
constexpr Scalar<T> norm(const Vector<T>& v) {
    return Scalar<T>(std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z));
}

template<typename T>
/**
 * @brief Computes the euclidian norm (magnitude) of a bivector in 3D space.
 * @param b The bivector for which to compute the norm.
 * @return The euclidian norm (magnitude) of the bivector b.
 */
constexpr Scalar<T> norm(const Bivector<T>& b) {
    return Scalar<T>(std::sqrt(b.xy * b.xy + b.xz * b.xz + b.yz * b.yz));  
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
 * @brief Computes the squared euclidian norm (magnitude) of a 3D vector.
 * @param v The vector for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the vector v.
 */
constexpr Scalar<T> squared_norm(const Vector<T>& v) {
    return Scalar<T>(v.x * v.x + v.y * v.y + v.z * v.z);
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a bivector in 3D space.
 * @param b The bivector for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the bivector b.
 */
constexpr Scalar<T> squared_norm(const Bivector<T>& b) {
    return Scalar<T>(b.xy * b.xy + b.xz * b.xz + b.yz * b.yz);
}

template<typename T>
/**
 * @brief Computes the euclidian norm (magnitude) of a trivector in 3D space.
 * @param t The trivector for which to compute the norm.
 * @return The euclidian norm (magnitude) of the trivector t.
 */
constexpr Scalar<T> norm(const Trivector<T>& t) {
    return Scalar<T>(std::abs(t.e123));
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a trivector in 3D space.
 * @param t The trivector for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the trivector t.
 */
constexpr Scalar<T> squared_norm(const Trivector<T>& t) {
    return Scalar<T>(t.e123 * t.e123);
}

template<typename T>
/**
 * @brief Computes the euclidian norm (magnitude) of a multivector in 3D space.
 * @param m The multivector for which to compute the norm.
 * @return The euclidian norm (magnitude) of the multivector m.
 */
constexpr Scalar<T> norm(const Multivector<T>& m) {
    return Scalar<T>(std::sqrt(squared_norm(m).value));
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a multivector in 3D space.
 * @param m The multivector for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the multivector m.
 */
constexpr Scalar<T> squared_norm(const Multivector<T>& m) {
    return Scalar<T>((m.scalar.value * m.scalar.value)
        + (m.vector.x * m.vector.x) + (m.vector.y * m.vector.y) + (m.vector.z * m.vector.z)
        + (m.bivector.xy * m.bivector.xy) + (m.bivector.xz * m.bivector.xz) + (m.bivector.yz * m.bivector.yz)
        + (m.trivector.e123 * m.trivector.e123));
}

template<typename T>
/**
 * @brief Computes the euclidian norm (magnitude) of a rotor in 3D space.
 * @param r The rotor for which to compute the norm.
 * @return The euclidian norm (magnitude) of the rotor r.
 */ 
constexpr Scalar<T> norm(const Rotor<T>& r) {
    return Scalar<T>(std::sqrt(squared_norm(r).value));
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a rotor in 3D space.
 * @param r The rotor for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the rotor r.
 */ 
constexpr Scalar<T> squared_norm(const Rotor<T>& r) {
    return Scalar<T>((r.scalar.value * r.scalar.value)
        + (r.bivector.xy * r.bivector.xy) + (r.bivector.xz * r.bivector.xz) + (r.bivector.yz * r.bivector.yz));
}

} // namespace CliffordCore::Cl3