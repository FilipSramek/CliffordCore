#pragma once

/**
 * @file norm.hpp
 * @brief Magnitudes: norm and squared_norm for every type.
 */

#include <cmath>
#include "../scalar.hpp"
#include "../vector3.hpp"
#include "../bivector3.hpp"
#include "../trivector3.hpp"
#include "../multivector3.hpp"
#include "../rotor3.hpp"

namespace CliffordCore
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
constexpr Scalar<T> norm(const Vector3<T>& v) {
    return Scalar<T>(std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z));
}

template<typename T>
/**
 * @brief Computes the euclidian norm (magnitude) of a bivector in 3D space.
 * @param b The bivector for which to compute the norm.
 * @return The euclidian norm (magnitude) of the bivector b.
 */
constexpr Scalar<T> norm(const Bivector3<T>& b) {
    return Scalar<T>(std::sqrt(b.xy * b.xy + b.xz * b.xz + b.yz * b.yz));  
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a 3D vector.
 * @param v The vector for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the vector v.
 */
constexpr Scalar<T> squared_norm(const Vector3<T>& v) {
    return Scalar<T>(v.x * v.x + v.y * v.y + v.z * v.z);
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a bivector in 3D space.
 * @param b The bivector for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the bivector b.
 */
constexpr Scalar<T> squared_norm(const Bivector3<T>& b) {
    return Scalar<T>(b.xy * b.xy + b.xz * b.xz + b.yz * b.yz);
}

template<typename T>
/**
 * @brief Computes the euclidian norm (magnitude) of a trivector in 3D space.
 * @param t The trivector for which to compute the norm.
 * @return The euclidian norm (magnitude) of the trivector t.
 */
constexpr Scalar<T> norm(const Trivector3<T>& t) {
    return Scalar<T>(std::abs(t.e123));
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a trivector in 3D space.
 * @param t The trivector for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the trivector t.
 */
constexpr Scalar<T> squared_norm(const Trivector3<T>& t) {
    return Scalar<T>(t.e123 * t.e123);
}

template<typename T>
/**
 * @brief Computes the euclidian norm (magnitude) of a multivector in 3D space.
 * @param m The multivector for which to compute the norm.
 * @return The euclidian norm (magnitude) of the multivector m.
 */
constexpr Scalar<T> norm(const Multivector3<T>& m) {
    return Scalar<T>(std::sqrt(squared_norm(m).value));
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a multivector in 3D space.
 * @param m The multivector for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the multivector m.
 */
constexpr Scalar<T> squared_norm(const Multivector3<T>& m) {
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
constexpr Scalar<T> norm(const Rotor3<T>& r) {
    return Scalar<T>(std::sqrt(squared_norm(r).value));
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a rotor in 3D space.
 * @param r The rotor for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the rotor r.
 */ 
constexpr Scalar<T> squared_norm(const Rotor3<T>& r) {
    return Scalar<T>((r.scalar.value * r.scalar.value)
        + (r.bivector.xy * r.bivector.xy) + (r.bivector.xz * r.bivector.xz) + (r.bivector.yz * r.bivector.yz));
}

} // namespace CliffordCore