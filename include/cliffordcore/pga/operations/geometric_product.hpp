#pragma once

/**
 * @file cliffordcore/pga/operations/geometric_product.hpp
 * @brief Geometric products, including the full 16x16 multiplication table.
 *
 * The full 16x16 Cayley table of Cl(3,0,1) lives here, as
 * geometric_product(Multivector, Multivector). Every mixed product in the
 * library routes through it, so there is exactly one multiplication table
 * to get right. The even-subalgebra closures -- Rotor * Rotor,
 * Translator * Translator and Motor * Motor -- are written out directly for
 * efficiency; the tests assert they agree with the table.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <type_traits>
#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../quadvector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"
#include "../translator.hpp"
#include "../motor.hpp"
#include "dot_product.hpp"
#include "wedge_product.hpp"

// Geometric products.
//
// Basis convention: e0 e1 e2 e3 with e0^2 = 0 and e1^2 = e2^2 = e3^2 = +1,
// every grade stored in lexicographic order: (e01, e02, e03, e12, e13, e23),
// (e012, e013, e023, e123), e0123. Of the 256 entries in the table, 64 vanish
// because the product contains e0 twice. The sign of each surviving entry is
// the parity of the permutation that sorts the concatenated indices, and the
// table below was generated from exactly that rule and checked against it in
// tests/test_pga.cpp -- if you edit a term by hand, the test will tell you.
//
// operator* always returns the most general type the product can produce, so
// there is exactly one operator* per operand pair. Named functions provide the
// narrower spellings -- motor_product() packs a plane-plane product into a
// Motor -- because overloads cannot differ by return type alone.

namespace CliffordCore::PGA
{
    template<typename T>
    /**
     * @brief Computes the geometric product of two planes, resulting in a multivector.
     * @param a The first vector.
     * @param b The second vector.
     * @return The resulting multivector: the dot product in grade 0 and the
     *         wedge product (their meet) in grade 2. Grades 1, 3 and 4 are zero.
     */
    constexpr Multivector<T> geometric_product(const Vector<T>& a, const Vector<T>& b) {
        return Multivector<T>(
            a | b,
            Vector<T>(),
            a ^ b,
            Trivector<T>(),
            Quadvector<T>()
        );
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two planes, resulting in a multivector.
     * @param a The first vector.
     * @param b The second vector.
     * @return The resulting multivector = dot_product(a, b) + wedge_product(a, b).
     */
    constexpr Multivector<T> operator*(const Vector<T>& a, const Vector<T>& b) {
        return geometric_product(a, b);
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two planes, packed as a motor.
     * @param a The first vector.
     * @param b The second vector.
     * @return The resulting motor: the composition of the reflections in b and
     *         then a, which is a rotation about their common line by twice the
     *         angle between them, or a translation if they are parallel. Same
     *         components as geometric_product(a, b), whose odd grades are zero.
     */
    constexpr Motor<T> motor_product(const Vector<T>& a, const Vector<T>& b) {
        return Motor<T>(a | b, a ^ b, Quadvector<T>());
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
        // Left operand components. `a` is the scalar; a0..a3 are the vector.
        const T a     = m.scalar.value;
        const T a0    = m.vector.e0,     a1    = m.vector.e1,     a2    = m.vector.e2,   a3 = m.vector.e3;
        const T a01   = m.bivector.e01,  a02   = m.bivector.e02,  a03   = m.bivector.e03;
        const T a12   = m.bivector.e12,  a13   = m.bivector.e13,  a23   = m.bivector.e23;
        const T a012  = m.trivector.e012, a013 = m.trivector.e013, a023 = m.trivector.e023, a123 = m.trivector.e123;
        const T a0123 = m.quadvector.e0123;

        // Right operand components.
        const T b     = n.scalar.value;
        const T b0    = n.vector.e0,     b1    = n.vector.e1,     b2    = n.vector.e2,   b3 = n.vector.e3;
        const T b01   = n.bivector.e01,  b02   = n.bivector.e02,  b03   = n.bivector.e03;
        const T b12   = n.bivector.e12,  b13   = n.bivector.e13,  b23   = n.bivector.e23;
        const T b012  = n.trivector.e012, b013 = n.trivector.e013, b023 = n.trivector.e023, b123 = n.trivector.e123;
        const T b0123 = n.quadvector.e0123;

        // Each coefficient is grouped by the grade of the left factor, one line
        // per grade. The non-e0 blades collect 8 terms, the e0 blades 16.

        // Grade 0. Only the Euclidean blades square to something non-zero.
        const T c = a*b
                    + a1*b1 + a2*b2 + a3*b3
                    - a12*b12 - a13*b13 - a23*b23
                    - a123*b123;   // 8 terms

        // Grade 1.
        const T c0 = a*b0
                     + a0*b - a1*b01 - a2*b02 - a3*b03
                     + a01*b1 + a02*b2 + a03*b3 - a12*b012 - a13*b013 - a23*b023
                     - a012*b12 - a013*b13 - a023*b23 + a123*b0123
                     - a0123*b123;   // 16 terms
        const T c1 = a*b1
                     + a1*b - a2*b12 - a3*b13
                     + a12*b2 + a13*b3 - a23*b123
                     - a123*b23;   // 8 terms
        const T c2 = a*b2
                     + a1*b12 + a2*b - a3*b23
                     - a12*b1 + a13*b123 + a23*b3
                     + a123*b13;   // 8 terms
        const T c3 = a*b3
                     + a1*b13 + a2*b23 + a3*b
                     - a12*b123 - a13*b1 - a23*b2
                     - a123*b12;   // 8 terms

        // Grade 2.
        const T c01 = a*b01
                      + a0*b1 - a1*b0 + a2*b012 + a3*b013
                      + a01*b - a02*b12 - a03*b13 + a12*b02 + a13*b03 - a23*b0123
                      + a012*b2 + a013*b3 - a023*b123 + a123*b023
                      - a0123*b23;   // 16 terms
        const T c02 = a*b02
                      + a0*b2 - a1*b012 - a2*b0 + a3*b023
                      + a01*b12 + a02*b - a03*b23 - a12*b01 + a13*b0123 + a23*b03
                      - a012*b1 + a013*b123 + a023*b3 - a123*b013
                      + a0123*b13;   // 16 terms
        const T c03 = a*b03
                      + a0*b3 - a1*b013 - a2*b023 - a3*b0
                      + a01*b13 + a02*b23 + a03*b - a12*b0123 - a13*b01 - a23*b02
                      - a012*b123 - a013*b1 - a023*b2 + a123*b012
                      - a0123*b12;   // 16 terms
        const T c12 = a*b12
                      + a1*b2 - a2*b1 + a3*b123
                      + a12*b - a13*b23 + a23*b13
                      + a123*b3;   // 8 terms
        const T c13 = a*b13
                      + a1*b3 - a2*b123 - a3*b1
                      + a12*b23 + a13*b - a23*b12
                      - a123*b2;   // 8 terms
        const T c23 = a*b23
                      + a1*b123 + a2*b3 - a3*b2
                      - a12*b13 + a13*b12 + a23*b
                      + a123*b1;   // 8 terms

        // Grade 3.
        const T c012 = a*b012
                       + a0*b12 - a1*b02 + a2*b01 - a3*b0123
                       + a01*b2 - a02*b1 + a03*b123 + a12*b0 - a13*b023 + a23*b013
                       + a012*b - a013*b23 + a023*b13 - a123*b03
                       + a0123*b3;   // 16 terms
        const T c013 = a*b013
                       + a0*b13 - a1*b03 + a2*b0123 + a3*b01
                       + a01*b3 - a02*b123 - a03*b1 + a12*b023 + a13*b0 - a23*b012
                       + a012*b23 + a013*b - a023*b12 + a123*b02
                       - a0123*b2;   // 16 terms
        const T c023 = a*b023
                       + a0*b23 - a1*b0123 - a2*b03 + a3*b02
                       + a01*b123 + a02*b3 - a03*b2 - a12*b013 + a13*b012 + a23*b0
                       - a012*b13 + a013*b12 + a023*b - a123*b01
                       + a0123*b1;   // 16 terms
        const T c123 = a*b123
                       + a1*b23 - a2*b13 + a3*b12
                       + a12*b3 - a13*b2 + a23*b1
                       + a123*b;   // 8 terms

        // Grade 4.
        const T c0123 = a*b0123
                        + a0*b123 - a1*b023 + a2*b013 - a3*b012
                        + a01*b23 - a02*b13 + a03*b12 + a12*b03 - a13*b02 + a23*b01
                        + a012*b3 - a013*b2 + a023*b1 - a123*b0
                        + a0123*b;   // 16 terms

        return Multivector<T>(
            Scalar<T>(c),
            Vector<T>(c0, c1, c2, c3),
            Bivector<T>(c01, c02, c03, c12, c13, c23),
            Trivector<T>(c012, c013, c023, c123),
            Quadvector<T>(c0123)
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
     * @brief Computes the geometric product of two rotors, resulting in a rotor. Analogous to two consecutive rotations.
     * @param r The first rotor.
     * @param s The second rotor.
     * @return The resulting rotor, equivalent to applying s then r. Rotors
     *         contain no e0, so this is the same closed form as in Cl3.
     */
    constexpr Rotor<T> operator*(const Rotor<T>& r, const Rotor<T>& s) {
        const T a = r.scalar.value, a12 = r.e12, a13 = r.e13, a23 = r.e23;
        const T b = s.scalar.value, b12 = s.e12, b13 = s.e13, b23 = s.e23;

        return Rotor<T>(
            Scalar<T>(a*b - a12*b12 - a13*b13 - a23*b23),
            a*b12 + a12*b - a13*b23 + a23*b13,
            a*b13 + a12*b23 + a13*b - a23*b12,
            a*b23 - a12*b13 + a13*b12 + a23*b
        );
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two rotors.
     * @param r The first rotor.
     * @param s The second rotor.
     * @return The resulting rotor, equivalent to applying s then r.
     */
    constexpr Rotor<T> rotor_product(const Rotor<T>& r, const Rotor<T>& s) {
        return r * s;
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two translators, resulting in a translator.
     * @param t The first translator.
     * @param u The second translator.
     * @return The resulting translator. Every ideal bivector squares to zero
     *         and any two of them multiply to zero, so
     *         (s1 + a)(s2 + b) = s1 s2 + s1 b + s2 a: for normalised
     *         translators the displacements simply add.
     */
    constexpr Translator<T> operator*(const Translator<T>& t, const Translator<T>& u) {
        const T a = t.scalar.value;
        const T b = u.scalar.value;
        return Translator<T>(
            Scalar<T>(a * b),
            a * u.e01 + b * t.e01,
            a * u.e02 + b * t.e02,
            a * u.e03 + b * t.e03
        );
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two translators.
     * @param t The first translator.
     * @param u The second translator.
     * @return The resulting translator, equivalent to applying u then t.
     */
    constexpr Translator<T> translator_product(const Translator<T>& t, const Translator<T>& u) {
        return t * u;
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two motors, resulting in a motor. Analogous to two consecutive rigid motions.
     * @param m The first motor.
     * @param n The second motor.
     * @return The resulting motor, equivalent to applying n then m. This is the
     *         even-even block of the 16x16 table above.
     */
    constexpr Motor<T> operator*(const Motor<T>& m, const Motor<T>& n) {
        const T a     = m.scalar.value;
        const T a01   = m.bivector.e01, a02 = m.bivector.e02, a03 = m.bivector.e03;
        const T a12   = m.bivector.e12, a13 = m.bivector.e13, a23 = m.bivector.e23;
        const T a0123 = m.quadvector.e0123;

        const T b     = n.scalar.value;
        const T b01   = n.bivector.e01, b02 = n.bivector.e02, b03 = n.bivector.e03;
        const T b12   = n.bivector.e12, b13 = n.bivector.e13, b23 = n.bivector.e23;
        const T b0123 = n.quadvector.e0123;

        const T c = a*b
                    - a12*b12 - a13*b13 - a23*b23;   // 4 terms
        const T c01 = a*b01
                      + a01*b - a02*b12 - a03*b13 + a12*b02 + a13*b03 - a23*b0123
                      - a0123*b23;   // 8 terms
        const T c02 = a*b02
                      + a01*b12 + a02*b - a03*b23 - a12*b01 + a13*b0123 + a23*b03
                      + a0123*b13;   // 8 terms
        const T c03 = a*b03
                      + a01*b13 + a02*b23 + a03*b - a12*b0123 - a13*b01 - a23*b02
                      - a0123*b12;   // 8 terms
        const T c12 = a*b12
                      + a12*b - a13*b23 + a23*b13;   // 4 terms
        const T c13 = a*b13
                      + a12*b23 + a13*b - a23*b12;   // 4 terms
        const T c23 = a*b23
                      - a12*b13 + a13*b12 + a23*b;   // 4 terms
        const T c0123 = a*b0123
                        + a01*b23 - a02*b13 + a03*b12 + a12*b03 - a13*b02 + a23*b01
                        + a0123*b;   // 8 terms

        return Motor<T>(
            Scalar<T>(c),
            Bivector<T>(c01, c02, c03, c12, c13, c23),
            Quadvector<T>(c0123)
        );
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two motors.
     * @param m The first motor.
     * @param n The second motor.
     * @return The resulting motor, equivalent to applying n then m.
     */
    constexpr Motor<T> motor_product(const Motor<T>& m, const Motor<T>& n) {
        return m * n;
    }
} // namespace CliffordCore::PGA
