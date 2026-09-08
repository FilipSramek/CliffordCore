#pragma once

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
    return Scalar<T>(std::abs(t.value));
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a trivector in 3D space.
 * @param t The trivector for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the trivector t.
 */
constexpr Scalar<T> squared_norm(const Trivector3<T>& t) {
    return Scalar<T>(t.value * t.value);
}

template<typename T>
/**
 * @brief Computes the euclidian norm (magnitude) of a multivector in 3D space.
 * @param m The multivector for which to compute the norm.
 * @return The euclidian norm (magnitude) of the multivector m.
 */
constexpr Scalar<T> norm(const Multivector3<T>& m) {
    return Scalar<T>(std::sqrt((m.scalar * m.scalar) + (m.vector.x * m.vector.x) + (m.vector.y * m.vector.y) + (m.vector.z * m.vector.z) + (m.bivector.xy * m.bivector.xy) + (m.bivector.xz * m.bivector.xz) + (m.bivector.yz * m.bivector.yz) + (m.trivector.value * m.trivector.value)));
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a multivector in 3D space.
 * @param m The multivector for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the multivector m.
 */
constexpr Scalar<T> squared_norm(const Multivector3<T>& m) {
    return Scalar<T>((m.scalar * m.scalar) + (m.vector.x * m.vector.x) + (m.vector.y * m.vector.y) + (m.vector.z * m.vector.z) + (m.bivector.xy * m.bivector.xy) + (m.bivector.xz * m.bivector.xz) + (m.bivector.yz * m.bivector.yz) + (m.trivector.value * m.trivector.value));
}

template<typename T>
/**
 * @brief Computes the euclidian norm (magnitude) of a rotor in 3D space.
 * @param r The rotor for which to compute the norm.
 * @return The euclidian norm (magnitude) of the rotor r.
 */ 
constexpr Scalar<T> norm(const Rotor3<T>& r) {
    return Scalar<T>(std::sqrt((r.scalar * r.scalar) + (r.bivector.xy * r.bivector.xy) + (r.bivector.xz * r.bivector.xz) + (r.bivector.yz * r.bivector.yz)));
}

template<typename T>
/**
 * @brief Computes the squared euclidian norm (magnitude) of a rotor in 3D space.
 * @param r The rotor for which to compute the squared norm.
 * @return The squared euclidian norm (magnitude) of the rotor r.
 */ 
constexpr Scalar<T> squared_norm(const Rotor3<T>& r) {
    return Scalar<T>((r.scalar * r.scalar) + (r.bivector.xy * r.bivector.xy) + (r.bivector.xz * r.bivector.xz) + (r.bivector.yz * r.bivector.yz));
}

} // namespace CliffordCore