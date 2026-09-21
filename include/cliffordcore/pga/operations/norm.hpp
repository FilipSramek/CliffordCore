#pragma once

/**
 * @file cliffordcore/pga/operations/norm.hpp
 * @brief Magnitudes: norm, squared_norm and ideal_norm for every type.
 *
 * In Cl(3,0,1) a norm is not the root of the sum of every squared
 * component, because e0 squares to zero: norm(x) is sqrt(<x ~x>_0), which
 * is the root of the sum of squares of the components that contain no e0.
 * That is the Euclidean size of an object -- the length of a plane's
 * normal, of a line's direction, a point's weight -- and it is zero for
 * every ideal element. The e0-containing components have their own size,
 * ideal_norm, which is the norm of the complement: a point's ideal norm is
 * its distance from the origin times its weight, a line's is the size of
 * its moment.
 *
 * squared_norm avoids the square root where the comparison does not need
 * it. For Bivector, Motor and Multivector the full product x ~x is a
 * "study number" a + b e0123 rather than a scalar; squared_norm is the a.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <cmath>
#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../quadvector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"
#include "../translator.hpp"
#include "../motor.hpp"

namespace CliffordCore::PGA
{

// -------------------------------------------------------------------------
// squared_norm: the sum of the squared non-e0 components
// -------------------------------------------------------------------------

template<typename T>
/**
 * @brief Computes the squared norm of a scalar.
 * @param s The scalar.
 * @return The scalar squared.
 */
constexpr Scalar<T> squared_norm(const Scalar<T>& s) {
    return Scalar<T>(s.value * s.value);
}

template<typename T>
/**
 * @brief Computes the squared norm of a plane: the squared length of its normal.
 * @param v The vector.
 * @return e1^2 + e2^2 + e3^2. The e0 component does not contribute.
 */
constexpr Scalar<T> squared_norm(const Vector<T>& v) {
    return Scalar<T>(v.e1 * v.e1 + v.e2 * v.e2 + v.e3 * v.e3);
}

template<typename T>
/**
 * @brief Computes the squared norm of a line: the squared length of its direction.
 * @param b The bivector.
 * @return e12^2 + e13^2 + e23^2. The ideal components do not contribute.
 */
constexpr Scalar<T> squared_norm(const Bivector<T>& b) {
    return Scalar<T>(b.e12 * b.e12 + b.e13 * b.e13 + b.e23 * b.e23);
}

template<typename T>
/**
 * @brief Computes the squared norm of a point: its squared weight.
 * @param t The trivector.
 * @return e123^2. The ideal components do not contribute.
 */
constexpr Scalar<T> squared_norm(const Trivector<T>& t) {
    return Scalar<T>(t.e123 * t.e123);
}

template<typename T>
/**
 * @brief Computes the squared norm of a quadvector.
 * @param q The quadvector.
 * @return Zero, always: e0123 contains e0.
 */
constexpr Scalar<T> squared_norm(const Quadvector<T>& q) {
    (void)q;
    return Scalar<T>(T(0));
}

template<typename T>
/**
 * @brief Computes the squared norm of a multivector.
 * @param m The multivector.
 * @return The sum of the squares of the eight non-e0 components, which is
 *         the scalar part of m * reverse(m).
 */
constexpr Scalar<T> squared_norm(const Multivector<T>& m) {
    return Scalar<T>(
        squared_norm(m.scalar).value
        + squared_norm(m.vector).value
        + squared_norm(m.bivector).value
        + squared_norm(m.trivector).value
    );
}

template<typename T>
/**
 * @brief Computes the squared norm of a rotor.
 * @param r The rotor.
 * @return scalar^2 + e12^2 + e13^2 + e23^2, which is r * reverse(r).
 */
constexpr Scalar<T> squared_norm(const Rotor<T>& r) {
    return Scalar<T>(r.scalar.value * r.scalar.value + r.e12 * r.e12 + r.e13 * r.e13 + r.e23 * r.e23);
}

template<typename T>
/**
 * @brief Computes the squared norm of a translator.
 * @param t The translator.
 * @return scalar^2, which is t * reverse(t): the ideal part contributes nothing.
 */
constexpr Scalar<T> squared_norm(const Translator<T>& t) {
    return Scalar<T>(t.scalar.value * t.scalar.value);
}

template<typename T>
/**
 * @brief Computes the squared norm of a motor.
 * @param m The motor.
 * @return scalar^2 + e12^2 + e13^2 + e23^2: the scalar part of the study
 *         number m * reverse(m).
 */
constexpr Scalar<T> squared_norm(const Motor<T>& m) {
    return Scalar<T>(squared_norm(m.scalar).value + squared_norm(m.bivector).value);
}

// -------------------------------------------------------------------------
// norm: the square root of squared_norm
// -------------------------------------------------------------------------

template<typename T>
/**
 * @brief Computes the norm (magnitude) of a scalar.
 * @param s The scalar.
 * @return The absolute value of s.
 */
constexpr Scalar<T> norm(const Scalar<T>& s) {
    return Scalar<T>(std::abs(s.value));
}

template<typename T>
/**
 * @brief Computes the norm of a plane: the length of its normal.
 * @param v The vector.
 * @return sqrt(e1^2 + e2^2 + e3^2). Zero for the ideal plane.
 */
constexpr Scalar<T> norm(const Vector<T>& v) {
    return Scalar<T>(std::sqrt(squared_norm(v).value));
}

template<typename T>
/**
 * @brief Computes the norm of a line: the length of its direction.
 * @param b The bivector.
 * @return sqrt(e12^2 + e13^2 + e23^2). Zero for an ideal line.
 */
constexpr Scalar<T> norm(const Bivector<T>& b) {
    return Scalar<T>(std::sqrt(squared_norm(b).value));
}

template<typename T>
/**
 * @brief Computes the norm of a point: its absolute weight.
 * @param t The trivector.
 * @return |e123|. Zero for an ideal point.
 */
constexpr Scalar<T> norm(const Trivector<T>& t) {
    return Scalar<T>(std::abs(t.e123));
}

template<typename T>
/**
 * @brief Computes the norm of a quadvector.
 * @param q The quadvector.
 * @return Zero, always. Use ideal_norm for its size.
 */
constexpr Scalar<T> norm(const Quadvector<T>& q) {
    (void)q;
    return Scalar<T>(T(0));
}

template<typename T>
/**
 * @brief Computes the norm of a multivector.
 * @param m The multivector.
 * @return sqrt(squared_norm(m)).
 */
constexpr Scalar<T> norm(const Multivector<T>& m) {
    return Scalar<T>(std::sqrt(squared_norm(m).value));
}

template<typename T>
/**
 * @brief Computes the norm of a rotor.
 * @param r The rotor.
 * @return sqrt(squared_norm(r)); 1 for a unit rotor.
 */
constexpr Scalar<T> norm(const Rotor<T>& r) {
    return Scalar<T>(std::sqrt(squared_norm(r).value));
}

template<typename T>
/**
 * @brief Computes the norm of a translator.
 * @param t The translator.
 * @return |scalar|; 1 for a normalised translator.
 */
constexpr Scalar<T> norm(const Translator<T>& t) {
    return Scalar<T>(std::abs(t.scalar.value));
}

template<typename T>
/**
 * @brief Computes the norm of a motor.
 * @param m The motor.
 * @return sqrt(squared_norm(m)); 1 for a unit motor. Note that a motor can
 *         have unit norm and still not satisfy m * reverse(m) = 1, because
 *         of the e0123 part of that product; see normalize.hpp.
 */
constexpr Scalar<T> norm(const Motor<T>& m) {
    return Scalar<T>(std::sqrt(squared_norm(m).value));
}

// -------------------------------------------------------------------------
// ideal_norm: the size of the e0-containing part
// -------------------------------------------------------------------------

template<typename T>
/**
 * @brief Computes the ideal norm of a scalar.
 * @param s The scalar.
 * @return Zero: a scalar has no e0 part.
 */
constexpr Scalar<T> ideal_norm(const Scalar<T>& s) {
    (void)s;
    return Scalar<T>(T(0));
}

template<typename T>
/**
 * @brief Computes the ideal norm of a plane.
 * @param v The vector.
 * @return |e0|: the plane's offset from the origin times the length of its normal.
 */
constexpr Scalar<T> ideal_norm(const Vector<T>& v) {
    return Scalar<T>(std::abs(v.e0));
}

template<typename T>
/**
 * @brief Computes the ideal norm of a line.
 * @param b The bivector.
 * @return sqrt(e01^2 + e02^2 + e03^2): the size of the line's moment.
 */
constexpr Scalar<T> ideal_norm(const Bivector<T>& b) {
    return Scalar<T>(std::sqrt(b.e01 * b.e01 + b.e02 * b.e02 + b.e03 * b.e03));
}

template<typename T>
/**
 * @brief Computes the ideal norm of a point.
 * @param t The trivector.
 * @return sqrt(e012^2 + e013^2 + e023^2): the point's distance from the
 *         origin times its weight.
 */
constexpr Scalar<T> ideal_norm(const Trivector<T>& t) {
    return Scalar<T>(std::sqrt(t.e012 * t.e012 + t.e013 * t.e013 + t.e023 * t.e023));
}

template<typename T>
/**
 * @brief Computes the ideal norm of a quadvector.
 * @param q The quadvector.
 * @return |e0123|.
 */
constexpr Scalar<T> ideal_norm(const Quadvector<T>& q) {
    return Scalar<T>(std::abs(q.e0123));
}

template<typename T>
/**
 * @brief Computes the ideal norm of a multivector.
 * @param m The multivector.
 * @return The root of the sum of the squares of the eight e0 components.
 */
constexpr Scalar<T> ideal_norm(const Multivector<T>& m) {
    return Scalar<T>(std::sqrt(
        m.vector.e0 * m.vector.e0
        + m.bivector.e01 * m.bivector.e01 + m.bivector.e02 * m.bivector.e02 + m.bivector.e03 * m.bivector.e03
        + m.trivector.e012 * m.trivector.e012 + m.trivector.e013 * m.trivector.e013 + m.trivector.e023 * m.trivector.e023
        + m.quadvector.e0123 * m.quadvector.e0123
    ));
}

template<typename T>
/**
 * @brief Computes the ideal norm of a rotor.
 * @param r The rotor.
 * @return Zero: a rotor has no e0 part.
 */
constexpr Scalar<T> ideal_norm(const Rotor<T>& r) {
    (void)r;
    return Scalar<T>(T(0));
}

template<typename T>
/**
 * @brief Computes the ideal norm of a translator.
 * @param t The translator.
 * @return sqrt(e01^2 + e02^2 + e03^2): half the displacement for a normalised translator.
 */
constexpr Scalar<T> ideal_norm(const Translator<T>& t) {
    return Scalar<T>(std::sqrt(t.e01 * t.e01 + t.e02 * t.e02 + t.e03 * t.e03));
}

template<typename T>
/**
 * @brief Computes the ideal norm of a motor.
 * @param m The motor.
 * @return The root of the sum of the squares of the ideal bivector components and e0123.
 */
constexpr Scalar<T> ideal_norm(const Motor<T>& m) {
    return Scalar<T>(std::sqrt(
        m.bivector.e01 * m.bivector.e01 + m.bivector.e02 * m.bivector.e02 + m.bivector.e03 * m.bivector.e03
        + m.quadvector.e0123 * m.quadvector.e0123
    ));
}
} // namespace CliffordCore::PGA
