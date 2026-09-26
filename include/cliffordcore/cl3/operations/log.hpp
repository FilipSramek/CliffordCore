#pragma once

/**
 * @file cliffordcore/cl3/operations/log.hpp
 * @brief Logarithm of a rotor, producing a bivector.
 *
 * The inverse of exp: recovers the bivector a rotor exponentiates from.
 * Used by slerp to interpolate in the plane of rotation.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <cmath>
#include "../rotor.hpp"
#include "../bivector.hpp"
#include "../scalar.hpp"
#include "../operations/norm.hpp"

namespace CliffordCore::Cl3
{
    template<typename T>
    /**
     * @brief Computes the logarithm of a rotor in 3D space, resulting in a bivector.
     * @param r The rotor for which to compute the logarithm.
     * @return The resulting bivector from the logarithm of the rotor r.
     */
    constexpr Bivector<T> log(const Rotor<T>& r) {
        const T magnitude = r.bivector.magnitude().value;

        if (magnitude == T(0)) { // Handle the case when the bivector is zero
            return Bivector<T>(0, 0, 0);
        }

        // The half-angle, from both parts of the rotor at once. This used to be
        // acos(scalar), which is badly conditioned near the identity: its slope
        // is 1/sin(angle), so a rotor turning by 1e-3 lost three digits, and
        // float lost the angle altogether. atan2 keeps full relative precision
        // at every angle, needs no clamp against rounding just outside [-1, 1],
        // and depends only on the ratio of the two parts, so a rotor that has
        // drifted off unit length still gives the right angle. Cl(2,0) and
        // Cl(3,0,1) compute their logs the same way.
        const T angle = std::atan2(magnitude, r.scalar.value);

        // The angle scaled onto the unit bivector that carries the rotation plane.
        const T scale = angle / magnitude;

        return Bivector<T>(
            scale * r.bivector.xy,
            scale * r.bivector.xz,
            scale * r.bivector.yz
        );
    }
} // namespace CliffordCore::Cl3
