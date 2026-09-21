#pragma once

/**
 * @file cliffordcore/pga/operations/exp.hpp
 * @brief Exponential of a bivector, producing a motor.
 *
 * A bivector in Cl(3,0,1) is a screw: an axis, an angle and a pitch. Its
 * exponential is the motor that performs that screw motion, and every
 * proper rigid motion of 3D space is exp of some bivector. Cl3's rotor
 * exp is the special case where the bivector is purely Euclidean.
 *
 * The series does not collapse to cos + sin as it does in Cl3, because a
 * general bivector B does not square to a scalar: B^2 = -theta^2 + 2 t I,
 * a "study number", where theta is the Euclidean norm and I = e0123. Since
 * I^2 = 0 the study numbers behave like dual numbers, and the closed form
 * carries a derivative term alongside each of cos and sin. Concretely,
 * with B_e the Euclidean part, B_i the ideal part and d the pitch:
 *
 *     theta = |B_e|,  d = -<B^2>_4 / (2 theta)
 *     exp(B) = cos(theta)
 *            + (sin(theta)/theta) B
 *            + (d/theta) (cos(theta) - sin(theta)/theta) (B_e I)
 *            - d sin(theta) I
 *
 * which reduces to Cl3's cos + sin B/|B| when B is Euclidean, and to the
 * exact 1 + B when B is ideal -- a translator, whose series terminates
 * after two terms because every ideal bivector squares to zero.
 *
 * The formula here was checked against a forty-term power series on random
 * bivectors before being written down; see docs/pga.md.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <cmath>
#include "../scalar.hpp"
#include "../bivector.hpp"
#include "../quadvector.hpp"
#include "../motor.hpp"
#include "../translator.hpp"
#include "norm.hpp"

namespace CliffordCore::PGA
{
    template<typename T>
    /**
     * @brief Computes the exponential of a bivector, resulting in a motor.
     * @param b The bivector (screw) for which to compute the exponential.
     * @return The motor performing that screw motion. Always a unit motor:
     *         m * reverse(m) = 1. A purely ideal bivector gives the translator
     *         1 + b exactly; a purely Euclidean one gives a rotation about a
     *         line through the origin, by twice the bivector's norm.
     */
    constexpr Motor<T> exp(const Bivector<T>& b) {
        const T theta = norm(b).value;

        if (theta == T(0)) {
            // Purely ideal: the series terminates, because B^2 = 0.
            return Motor<T>(Scalar<T>(1), b, Quadvector<T>());
        }

        // <B^2>_4 = 2 (e01 e23 - e02 e13 + e03 e12); d is the pitch, the
        // distance travelled along the axis per unit of rotation angle.
        const T b4 = T(2) * (b.e01 * b.e23 - b.e02 * b.e13 + b.e03 * b.e12);
        const T d = -b4 / (T(2) * theta);

        const T c = std::cos(theta);
        const T s = std::sin(theta) / theta;

        // B_e * e0123 maps (e12, e13, e23) onto (-e03, +e02, -e01).
        const T k = (d / theta) * (c - s);

        return Motor<T>(
            Scalar<T>(c),
            Bivector<T>(
                s * b.e01 - k * b.e23,
                s * b.e02 + k * b.e13,
                s * b.e03 - k * b.e12,
                s * b.e12,
                s * b.e13,
                s * b.e23
            ),
            Quadvector<T>(-d * std::sin(theta))
        );
    }
} // namespace CliffordCore::PGA
