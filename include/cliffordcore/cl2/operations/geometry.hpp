#pragma once

/**
 * @file cliffordcore/cl2/operations/geometry.hpp
 * @brief Reflection, projection and rejection, plus the planar extras:
 *        perp, signed area and angles.
 *
 * reflect, project and reject, written with the contraction and the inverse
 * rather than assuming unit-length arguments; the reflection is -n v n^-1.
 * Alongside them, the operations that are specific to working in a plane:
 * the perpendicular, the signed area two vectors span, and angles read back
 * with atan2.
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
#include "dot_product.hpp"
#include "dual.hpp"
#include "geometric_product.hpp"
#include "grade.hpp"
#include "inverse.hpp"
#include "mixed_products.hpp"
#include "wedge_product.hpp"

// The classical geometry answers, as named functions.
//
// They divide through by inverse() rather than assuming a unit argument, so
// there is no silent precondition: pass any non-zero vector. Projection and
// rejection split a vector into the part along a direction and the part
// across it, so project(v, u) + reject(v, u) == v always.
//
// Two things Cl(3,0) has are deliberately missing. There are no overloads
// taking a Bivector: in two dimensions the bivector is the whole plane, so
// projecting onto it is the identity, rejecting from it gives zero, and
// reflecting in it means nothing. And a line through the origin can be named
// either by its direction or by its normal, both of which are Vectors, so only
// one of them can be reflect(Vector, Vector). This one takes the NORMAL, to
// match Cl(3,0)'s "the plane perpendicular to n". Reflecting in the line ALONG
// u is -reflect(v, u), since the two mirrors are a quarter turn apart.
//
// The angle functions return a raw T rather than a Scalar, like Cl(3,0,1)'s
// distance() and angle(), and are inline rather than constexpr because they
// call std::atan2.
//
// They promise the principal range (-pi, pi], and plain atan2 does not keep
// that promise: atan2(-0.0, -1) is -pi, and -0.0 is exactly what negating a
// vector with a zero component produces -- angle(-e1) would come back as -pi.
// Adding +0.0 to the first argument turns -0.0 into +0.0 and leaves every
// other value alone (IEEE 754 fixes that; a compiler may not fold it away
// without -ffast-math), so a half turn always reports +pi. log() and
// rotor_angle() do the same.

namespace CliffordCore::Cl2
{
    // -----------------------------------------------------------------------
    // Reflection
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Reflects a vector in the line through the origin perpendicular to n.
     * @param v The vector to reflect.
     * @param n The mirror line's normal. Need not be unit length.
     * @return The reflected vector, -n v inverse(n).
     *
     * Equivalent to the familiar v - 2*(v.n)*n for unit n, but written as a
     * product so it composes: two reflections make a rotation, which is what a
     * rotor is. To reflect in the line along a direction u instead, negate:
     * -reflect(v, u).
     */
    constexpr Vector<T> reflect(const Vector<T>& v, const Vector<T>& n) {
        return grade1(-(n * v * inverse(n)));
    }

    // -----------------------------------------------------------------------
    // Projection and rejection
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Projects a vector onto the line spanned by another vector.
     * @param v The vector to project.
     * @param u The direction to project onto. Need not be unit length.
     * @return The component of v parallel to u.
     */
    constexpr Vector<T> project(const Vector<T>& v, const Vector<T>& u) {
        return (v | u) * inverse(u);
    }

    template<typename T>
    /**
     * @brief Removes from a vector the part parallel to another vector.
     * @param v The vector to reject.
     * @param u The direction to reject from. Need not be unit length.
     * @return The component of v perpendicular to u. A bivector times a vector
     *         is already a pure vector in two dimensions, so no grade
     *         projection is needed.
     */
    constexpr Vector<T> reject(const Vector<T>& v, const Vector<T>& u) {
        return (v ^ u) * inverse(u);
    }

    // -----------------------------------------------------------------------
    // Planar extras
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief The vector turned a quarter turn counter-clockwise.
     * @param v The vector to turn.
     * @return (-v.y, v.x). This is exactly dual(v): grade 1 is its own dual
     *         in two dimensions, and multiplying by e12 is the quarter turn.
     */
    constexpr Vector<T> perp(const Vector<T>& v) {
        return dual(v);
    }

    template<typename T>
    /**
     * @brief The signed area of the parallelogram two vectors span.
     * @param a The first edge.
     * @param b The second edge.
     * @return a.x * b.y - a.y * b.x: positive when b lies counter-clockwise of
     *         a. This is the component of a ^ b -- the 2D stand-in for the
     *         cross product -- handed back as a Scalar. Halve it for a triangle.
     */
    constexpr Scalar<T> signed_area(const Vector<T>& a, const Vector<T>& b) {
        return Scalar<T>((a ^ b).xy);
    }

    template<typename T>
    /**
     * @brief The signed angle that turns one vector onto another.
     * @param a The starting direction. Need not be unit length.
     * @param b The target direction. Need not be unit length.
     * @return The angle in radians, in (-pi, pi], counter-clockwise positive:
     *         angle_between(e1, e2) == +pi/2. Deliberately NOT named angle():
     *         Cl(3,0,1)'s two-argument angle() is unsigned, and a silently
     *         different sign convention under a namespace switch is worse than a
     *         different name.
     */
    inline T angle_between(const Vector<T>& a, const Vector<T>& b) {
        // + T(0) turns a -0.0 into +0.0; see the note at the top of the file.
        return std::atan2((a ^ b).xy + T(0), (a | b).value);
    }

    template<typename T>
    /**
     * @brief The direction of a vector, as an angle from the positive x axis.
     * @param v The vector. Need not be unit length.
     * @return atan2(v.y, v.x), in (-pi, pi], counter-clockwise positive.
     */
    inline T angle(const Vector<T>& v) {
        // + T(0) turns a -0.0 into +0.0, so angle(-e1) is +pi, not -pi.
        return std::atan2(v.y + T(0), v.x);
    }
} // namespace CliffordCore::Cl2
