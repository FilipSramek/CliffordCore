#pragma once

/**
 * @file geometry.hpp
 * @brief Reflection, projection and rejection.
 *
 * reflect, project and reject, written with the contraction and the inverse
 * rather than assuming unit-length arguments. Projection onto a plane is
 * (v _| B) * inverse(B); the reflection is -n v n^-1.
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
#include "../multivector3.hpp"
#include "contraction.hpp"
#include "dual.hpp"
#include "geometric_product.hpp"
#include "grade.hpp"
#include "inverse.hpp"
#include "mixed_products.hpp"
#include "wedge_product.hpp"

// The classical geometry answers, as named functions.
//
// Each of these was previously spelled out by hand in the examples. They divide
// through by inverse() rather than assuming a unit argument, so there is no
// silent precondition: pass any non-zero vector or plane.
//
// Projection and rejection split a vector into the part inside a subspace and
// the part outside it, so project(v, x) + reject(v, x) == v always.

namespace CliffordCore
{
    // -----------------------------------------------------------------------
    // Reflection
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Reflects a vector in the plane through the origin perpendicular to n.
     * @param v The vector to reflect.
     * @param n The plane's normal. Need not be unit length.
     * @return The reflected vector, -n v inverse(n).
     *
     * Equivalent to the familiar v - 2*(v.n)*n for unit n, but written as a
     * product so it composes: two reflections make a rotation, which is what a
     * rotor is.
     */
    constexpr Vector3<T> reflect(const Vector3<T>& v, const Vector3<T>& n) {
        return grade1(-(n * v * inverse(n)));
    }

    template<typename T>
    /**
     * @brief Reflects a vector in a plane given as a bivector.
     * @param v The vector to reflect.
     * @param plane The plane to reflect in. Need not be unit length.
     * @return The reflected vector.
     *
     * The bivector names the plane directly rather than by its normal, which is
     * the formulation that carries over to algebras where a plane has no single
     * normal vector.
     */
    constexpr Vector3<T> reflect(const Vector3<T>& v, const Bivector3<T>& plane) {
        return reflect(v, dual(plane));
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
    constexpr Vector3<T> project(const Vector3<T>& v, const Vector3<T>& u) {
        return (v | u) * inverse(u);
    }

    template<typename T>
    /**
     * @brief Projects a vector onto a plane.
     * @param v The vector to project.
     * @param plane The plane to project onto. Need not be unit length.
     * @return The component of v lying in the plane.
     *
     * This is the one that needs the left contraction: (v _| B) * inverse(B).
     */
    constexpr Vector3<T> project(const Vector3<T>& v, const Bivector3<T>& plane) {
        return grade1((v | plane) * inverse(plane));
    }

    template<typename T>
    /**
     * @brief Removes from a vector the part parallel to another vector.
     * @param v The vector to reject.
     * @param u The direction to reject from. Need not be unit length.
     * @return The component of v perpendicular to u.
     */
    constexpr Vector3<T> reject(const Vector3<T>& v, const Vector3<T>& u) {
        return grade1((v ^ u) * inverse(u));
    }

    template<typename T>
    /**
     * @brief Removes from a vector the part lying in a plane.
     * @param v The vector to reject.
     * @param plane The plane to reject from. Need not be unit length.
     * @return The component of v perpendicular to the plane.
     */
    constexpr Vector3<T> reject(const Vector3<T>& v, const Bivector3<T>& plane) {
        return grade1((v ^ plane) * inverse(plane));
    }
} // namespace CliffordCore
