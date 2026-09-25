#pragma once

/**
 * @file cliffordcore/cl2/operations/exp.hpp
 * @brief Exponential of a bivector, producing a rotor.
 *
 * exp(theta * e12) = cos(theta) + sin(theta) e12 -- Euler's formula, with
 * e12 in the role of i. As in Cl(3,0), exp(theta * B) is the rotor rotating
 * by 2 * theta, so exp(pi/4 * e12) is a quarter turn -- clockwise, for the
 * reason given below.
 *
 * @author Filip Sramek
 * @version 0.3.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <cmath>
#include "../rotor.hpp"
#include "../bivector.hpp"
#include "../scalar.hpp"

// There is only one plane in two dimensions, so a bivector is a single signed
// number times e12, and e12^2 = -1 makes the exponential series collapse to
// cos + sin exactly as it does for complex numbers. No normalisation, and no
// special case for the zero bivector: cos(0) = 1 and sin(0) = 0 already give
// the identity rotor.
//
// Cl(3,0)'s exp computes |B| and then sin(|B|) * B/|B|. For B = theta e12 that
// is sin(|theta|) * theta/|theta| = sin(theta), since sine is odd -- so the two
// algebras agree, and this one simply does not need the detour.
//
// Orientation: sandwich(v, exp(theta e12)) turns v by -2 theta, i.e.
// CLOCKWISE for positive theta. That is why rotor_from_angle(angle) is
// exp(-(angle/2) e12). The sign is the same as Cl(3,0)'s about +z.

namespace CliffordCore::Cl2
{
    template<typename T>
    /**
     * @brief Computes the exponential of a bivector in 2D space, resulting in a rotor.
     * @param b The bivector for which to compute the exponential.
     * @return The resulting rotor, cos(b.xy) + sin(b.xy) e12. Always unit length.
     */
    constexpr Rotor<T> exp(const Bivector<T>& b) {
        return Rotor<T>(
            Scalar<T>(std::cos(b.xy)),
            Bivector<T>(std::sin(b.xy))
        );
    }
} // namespace CliffordCore::Cl2
