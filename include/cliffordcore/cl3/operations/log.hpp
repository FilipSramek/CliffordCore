#pragma once

/**
 * @file log.hpp
 * @brief Logarithm of a rotor, producing a bivector.
 *
 * The inverse of exp: recovers the bivector a rotor exponentiates from.
 * Used by slerp to interpolate in the plane of rotation.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <cmath>
#include "../rotor3.hpp"
#include "../bivector3.hpp"
#include "../scalar.hpp"
#include "../operations/norm.hpp"

namespace CliffordCore
{
    template<typename T>
    /**
     * @brief Computes the logarithm of a rotor in 3D space, resulting in a bivector.
     * @param r The rotor for which to compute the logarithm.
     * @return The resulting bivector from the logarithm of the rotor r.
     */
    constexpr Bivector3<T> log(const Rotor3<T>& r) {
        // A unit rotor has |scalar| <= 1, but rounding in exp() or in repeated
        // composition can leave it a few ulps outside that range, and std::acos
        // returns NaN off-domain. Clamp before the call.
        const T clamped_scalar = r.scalar.value < T(-1) ? T(-1)
                               : (r.scalar.value > T(1) ? T(1) : r.scalar.value);

        Scalar<T> arccos_scalar = Scalar<T> (std::acos(clamped_scalar));

        Scalar<T> Magnitude_bivector = r.bivector.magnitude();

        if (Magnitude_bivector.value == 0) { // Handle the case when the bivector is zero
            return Bivector3<T>(0, 0, 0);
        }

        // The angle scaled onto the unit bivector that carries the rotation plane.
        const T scale = arccos_scalar.value / Magnitude_bivector.value;

        return Bivector3<T>(
            scale * r.bivector.xy,
            scale * r.bivector.xz,
            scale * r.bivector.yz
        );
    }
} // namespace CliffordCore
