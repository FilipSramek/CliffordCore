#pragma once

/**
 * @file geometric_product.hpp
 * @brief Geometric products, including the full 8x8 multiplication table.
 */

#include <type_traits>
#include "../vector3.hpp"
#include "../multivector3.hpp"
#include "../rotor3.hpp"
#include "../bivector3.hpp"
#include "../scalar.hpp"
#include "dot_product.hpp"
#include "wedge_product.hpp"

// Geometric products.
//
// Basis convention: e1 e2 e3 with e_i^2 = +1, bivectors stored as
// (xy, xz, yz) = (e1e2, e1e3, e2e3), pseudoscalar e123 = e1e2e3.
// From that: e12^2 = e13^2 = e23^2 = -1 and e123^2 = -1.
//
// operator* always returns the most general type the product can produce, so
// there is exactly one operator* per operand pair. Named functions provide the
// narrower spellings -- rotor_product() packs a vector-vector product into a
// Rotor3 -- because overloads cannot differ by return type alone.

namespace CliffordCore
{
    template<typename T>
    /**
     * @brief Computes the geometric product of two 3D vectors, resulting in a multivector.
     * @param a The first vector.
     * @param b The second vector.
     * @return The resulting multivector from the geometric product = dot_product(a, b) + wedge_product(a, b) so that the scalar part is the dot product and the bivector part is the wedge product.
     */
    constexpr Multivector3<T> geometric_product(const Vector3<T>& a, const Vector3<T>& b) {
        return Multivector3<T>(
            a | b,
            Vector3<T>(0, 0, 0),
            a ^ b,
            Trivector3<T>(0)
        );
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two 3D vectors, resulting in a multivector.
     * @param a The first vector.
     * @param b The second vector.
     * @return The resulting multivector from the geometric product = dot_product(a, b) + wedge_product(a, b) so that the scalar part is the dot product and the bivector part is the wedge product.
     */
    constexpr Multivector3<T> operator*(const Vector3<T>& a, const Vector3<T>& b) {
        return geometric_product(a, b);
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two 3D vectors, packed as a rotor.
     * @param a The first vector.
     * @param b The second vector.
     * @return The resulting rotor. Same components as geometric_product(a, b),
     *         whose vector and trivector parts are always zero.
     */
    constexpr Rotor3<T> rotor_product(const Vector3<T>& a, const Vector3<T>& b) {
        return Rotor3<T>(
            a | b,
            a ^ b
        );
    }

    template<typename T>
    /**
     * @brief Computes the full geometric product of two multivectors.
     * @param m The left multivector.
     * @param n The right multivector.
     * @return The resulting multivector. This is the general product of the
     *         algebra; every other product here is a special case of it.
     */
    constexpr Multivector3<T> geometric_product(const Multivector3<T>& m, const Multivector3<T>& n) {
        // Left operand components.
        const T a0   = m.scalar.value;
        const T a1   = m.vector.x,    a2   = m.vector.y,  a3   = m.vector.z;
        const T a12  = m.bivector.xy, a13  = m.bivector.xz, a23 = m.bivector.yz;
        const T a123 = m.trivector.e123;

        // Right operand components.
        const T b0   = n.scalar.value;
        const T b1   = n.vector.x,    b2   = n.vector.y,  b3   = n.vector.z;
        const T b12  = n.bivector.xy, b13  = n.bivector.xz, b23 = n.bivector.yz;
        const T b123 = n.trivector.e123;

        // Grade 0. Bivectors and the pseudoscalar square to -1, hence the signs.
        const T c0 = a0*b0 + a1*b1 + a2*b2 + a3*b3
                   - a12*b12 - a13*b13 - a23*b23 - a123*b123;

        // Grade 1.
        const T c1 = a0*b1 + a1*b0
                   - a2*b12 - a3*b13 + a12*b2 + a13*b3
                   - a23*b123 - a123*b23;
        const T c2 = a0*b2 + a2*b0
                   + a1*b12 - a3*b23 - a12*b1 + a23*b3
                   + a13*b123 + a123*b13;
        const T c3 = a0*b3 + a3*b0
                   + a1*b13 + a2*b23 - a13*b1 - a23*b2
                   - a12*b123 - a123*b12;

        // Grade 2.
        const T c12 = a0*b12 + a12*b0
                    + a1*b2 - a2*b1
                    - a13*b23 + a23*b13
                    + a3*b123 + a123*b3;
        const T c13 = a0*b13 + a13*b0
                    + a1*b3 - a3*b1
                    + a12*b23 - a23*b12
                    - a2*b123 - a123*b2;
        const T c23 = a0*b23 + a23*b0
                    + a2*b3 - a3*b2
                    - a12*b13 + a13*b12
                    + a1*b123 + a123*b1;

        // Grade 3.
        const T c123 = a0*b123 + a123*b0
                     + a1*b23 - a2*b13 + a3*b12
                     + a23*b1 - a13*b2 + a12*b3;

        return Multivector3<T>(
            Scalar<T>(c0),
            Vector3<T>(c1, c2, c3),
            Bivector3<T>(c12, c13, c23),
            Trivector3<T>(c123)
        );
    }

    template<typename T>
    /**
     * @brief Computes the full geometric product of two multivectors.
     * @param m The left multivector.
     * @param n The right multivector.
     * @return The resulting multivector.
     */
    constexpr Multivector3<T> operator*(const Multivector3<T>& m, const Multivector3<T>& n) {
        return geometric_product(m, n);
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two 3D rotors, resulting in a 3D rotor. Analogous to two consecutive rotations.
     * @param r The first rotor.
     * @param s The second rotor.
     * @return The resulting rotor from the geometric product.
     */
    constexpr Rotor3<T> operator*(const Rotor3<T>& r, const Rotor3<T>& s) {
        const T r0 = r.scalar.value;
        const T r12 = r.bivector.xy, r13 = r.bivector.xz, r23 = r.bivector.yz;
        const T s0 = s.scalar.value;
        const T s12 = s.bivector.xy, s13 = s.bivector.xz, s23 = s.bivector.yz;

        return Rotor3<T>(
            Scalar<T>(r0*s0 - r12*s12 - r13*s13 - r23*s23),
            Bivector3<T>(
                r0*s12 + r12*s0 - r13*s23 + r23*s13,
                r0*s13 + r13*s0 + r12*s23 - r23*s12,
                r0*s23 + r23*s0 - r12*s13 + r13*s12
            )
        );
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two 3D rotors.
     * @param r The first rotor.
     * @param s The second rotor.
     * @return The resulting rotor, equivalent to applying s then r.
     */
    constexpr Rotor3<T> rotor_product(const Rotor3<T>& r, const Rotor3<T>& s) {
        return r * s;
    }
} // namespace CliffordCore
