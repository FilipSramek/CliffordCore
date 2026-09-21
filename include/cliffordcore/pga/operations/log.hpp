#pragma once

/**
 * @file cliffordcore/pga/operations/log.hpp
 * @brief Logarithm of a motor, producing a bivector.
 *
 * The inverse of exp: recovers the screw a motor performs, as a bivector
 * whose Euclidean norm is half the rotation angle and whose axis is the
 * screw axis. Used by slerp to interpolate motions, and by the caller who
 * wants the axis and angle back out of a motor.
 *
 * With s the motor's scalar part, B its bivector part, p its e0123 part
 * and n = |B_e|:
 *
 *     theta = atan2(n, s),  d = -p / n
 *     log(m) = (theta/n) B + (d/n) (1 - theta cos(theta)/n) (B_e I)
 *
 * A motor with no Euclidean bivector part is a translator, whose log is
 * just B / s. As in Cl3, the round trip is exact for angles inside the
 * principal range and wraps outside it.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <cmath>
#include "../scalar.hpp"
#include "../bivector.hpp"
#include "../motor.hpp"
#include "../rotor.hpp"
#include "../translator.hpp"
#include "norm.hpp"

namespace CliffordCore::PGA
{
    template<typename T>
    /**
     * @brief Computes the logarithm of a motor, resulting in a bivector.
     * @param m The motor for which to compute the logarithm. Expected to be a
     *        unit motor; normalize() it first if it has drifted.
     * @return The bivector (screw) with exp of it equal to m. For a pure
     *         translator this is the ideal bivector b with m = 1 + b.
     */
    constexpr Bivector<T> log(const Motor<T>& m) {
        const Bivector<T>& b = m.bivector;
        const T n = norm(b).value;

        if (n == T(0)) {
            // No rotation: a translator, whose log is the ideal part scaled
            // back by the scalar. A zero scalar too would mean a zero motor,
            // which has no logarithm; return the ideal part unchanged rather
            // than dividing by zero.
            const T s = m.scalar.value;
            if (s == T(0)) {
                return b;
            }
            return b / Scalar<T>(s);
        }

        // A unit motor has |scalar| <= 1, but rounding in exp() or in repeated
        // composition can leave it a few ulps outside; atan2 is well behaved
        // there, unlike acos, so no clamp is needed.
        const T theta = std::atan2(n, m.scalar.value);
        const T d = -m.quadvector.e0123 / n;

        // Inverse of the exp() coupling: B_e * e0123 maps
        // (e12, e13, e23) onto (-e03, +e02, -e01).
        const T k = (d / n) * (T(1) - theta * std::cos(theta) / n);
        const T scale = theta / n;

        return Bivector<T>(
            scale * b.e01 - k * b.e23,
            scale * b.e02 + k * b.e13,
            scale * b.e03 - k * b.e12,
            scale * b.e12,
            scale * b.e13,
            scale * b.e23
        );
    }

    template<typename T>
    /**
     * @brief Computes the logarithm of a rotor, resulting in a bivector.
     * @param r The rotor for which to compute the logarithm.
     * @return The Euclidean bivector with exp of it equal to r. This is the
     *         Cl3 case: no ideal part in, none out.
     */
    constexpr Bivector<T> log(const Rotor<T>& r) {
        return log(Motor<T>(r));
    }

    template<typename T>
    /**
     * @brief Computes the logarithm of a translator, resulting in a bivector.
     * @param t The translator for which to compute the logarithm.
     * @return The ideal bivector with exp of it equal to t, which is t's ideal
     *         part divided by its scalar.
     */
    constexpr Bivector<T> log(const Translator<T>& t) {
        return log(Motor<T>(t));
    }
} // namespace CliffordCore::PGA
