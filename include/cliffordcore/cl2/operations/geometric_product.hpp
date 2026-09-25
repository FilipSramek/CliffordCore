#pragma once

/**
 * @file cliffordcore/cl2/operations/geometric_product.hpp
 * @brief The geometric product: the full 4x4 table, plus fast paths.
 *
 * geometric_product(Multivector, Multivector) is the one multiplication
 * table in Cl(2,0); every mixed-grade product routes through it. Vector
 * times vector and rotor times rotor keep direct implementations for speed,
 * and the tests assert they agree with the general product.
 *
 * @author Filip Sramek
 * @version 0.3.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <type_traits>
#include "../vector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"
#include "../bivector.hpp"
#include "../scalar.hpp"
#include "dot_product.hpp"
#include "wedge_product.hpp"

// Geometric products.
//
// Basis convention: e1 e2 with e_i^2 = +1, and the one bivector
// e12 = e1e2, which is also the pseudoscalar. From that, e12^2 = -1.
//
// The whole table, verified against the permutation rule by
// test_cayley_table_against_permutation_rule:
//
//          1      e1     e2     e12
//   1   |  1      e1     e2     e12
//   e1  |  e1     1      e12    e2
//   e2  |  e2    -e12    1     -e1
//   e12 |  e12   -e2     e1    -1
//
// Note e1 e12 = +e2 but e12 e1 = -e2: the pseudoscalar ANTICOMMUTES with
// vectors here, where Cl(3,0)'s e123 commutes with everything. Most of what
// is different about this algebra follows from that one line.
//
// operator* always returns the most general type the product can produce, so
// there is exactly one operator* per operand pair. Named functions provide the
// narrower spellings -- rotor_product() packs a vector-vector product into a
// Rotor -- because overloads cannot differ by return type alone.

namespace CliffordCore::Cl2
{
    template<typename T>
    /**
     * @brief Computes the geometric product of two 2D vectors, resulting in a multivector.
     * @param a The first vector.
     * @param b The second vector.
     * @return The resulting multivector from the geometric product = dot_product(a, b) + wedge_product(a, b) so that the scalar part is the dot product and the bivector part is the wedge product.
     */
    constexpr Multivector<T> geometric_product(const Vector<T>& a, const Vector<T>& b) {
        return Multivector<T>(
            a | b,
            Vector<T>(0, 0),
            a ^ b
        );
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two 2D vectors, resulting in a multivector.
     * @param a The first vector.
     * @param b The second vector.
     * @return The resulting multivector. Its vector part is always zero, so the
     *         result is exactly a rotor -- but it is returned as a Multivector
     *         to match Cl(3,0). Use rotor_product(a, b) for the Rotor spelling.
     */
    constexpr Multivector<T> operator*(const Vector<T>& a, const Vector<T>& b) {
        return geometric_product(a, b);
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two 2D vectors, packed as a rotor.
     * @param a The first vector.
     * @param b The second vector.
     * @return The resulting rotor. Same components as geometric_product(a, b),
     *         whose vector part is always zero.
     */
    constexpr Rotor<T> rotor_product(const Vector<T>& a, const Vector<T>& b) {
        return Rotor<T>(
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
    constexpr Multivector<T> geometric_product(const Multivector<T>& m, const Multivector<T>& n) {
        // Left operand components.
        const T a0  = m.scalar.value;
        const T a1  = m.vector.x, a2 = m.vector.y;
        const T a12 = m.bivector.xy;

        // Right operand components.
        const T b0  = n.scalar.value;
        const T b1  = n.vector.x, b2 = n.vector.y;
        const T b12 = n.bivector.xy;

        // Four terms per coefficient, one from each row of the table above.

        // Grade 0. The bivector squares to -1, hence the sign.
        const T c0 = a0*b0 + a1*b1 + a2*b2 - a12*b12;

        // Grade 1. The e12 terms carry the anticommutation: e2 e12 = -e1 but
        // e12 e2 = +e1, and e1 e12 = +e2 but e12 e1 = -e2.
        const T c1 = a0*b1 + a1*b0 - a2*b12 + a12*b2;
        const T c2 = a0*b2 + a2*b0 + a1*b12 - a12*b1;

        // Grade 2.
        const T c12 = a0*b12 + a12*b0 + a1*b2 - a2*b1;

        return Multivector<T>(
            Scalar<T>(c0),
            Vector<T>(c1, c2),
            Bivector<T>(c12)
        );
    }

    template<typename T>
    /**
     * @brief Computes the full geometric product of two multivectors.
     * @param m The left multivector.
     * @param n The right multivector.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(const Multivector<T>& m, const Multivector<T>& n) {
        return geometric_product(m, n);
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two 2D rotors, resulting in a 2D rotor. Analogous to two consecutive rotations.
     * @param r The first rotor.
     * @param s The second rotor.
     * @return The resulting rotor from the geometric product.
     *
     * This is complex multiplication, term for term: (a + b i)(c + d i) =
     * (ac - bd) + (ad + bc) i with i = e12. The even subalgebra of Cl(2,0) is
     * commutative, so unlike in Cl(3,0), r * s == s * r here.
     */
    constexpr Rotor<T> operator*(const Rotor<T>& r, const Rotor<T>& s) {
        const T r0 = r.scalar.value, r12 = r.bivector.xy;
        const T s0 = s.scalar.value, s12 = s.bivector.xy;

        return Rotor<T>(
            Scalar<T>(r0*s0 - r12*s12),
            Bivector<T>(r0*s12 + r12*s0)
        );
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two 2D rotors.
     * @param r The first rotor.
     * @param s The second rotor.
     * @return The resulting rotor, equivalent to applying s then r.
     */
    constexpr Rotor<T> rotor_product(const Rotor<T>& r, const Rotor<T>& s) {
        return r * s;
    }
} // namespace CliffordCore::Cl2
