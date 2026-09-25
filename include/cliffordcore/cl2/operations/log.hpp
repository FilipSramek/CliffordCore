#pragma once

/**
 * @file cliffordcore/cl2/operations/log.hpp
 * @brief Logarithm of a rotor, producing a bivector.
 *
 * The inverse of exp: recovers the bivector a rotor exponentiates from.
 * Used by slerp to interpolate the angle of rotation. In Cl(2,0) this is the
 * argument of a complex number, computed with atan2.
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

// A rotor is a + b e12, a complex number, and its log is the angle it makes
// with the positive real axis: atan2(b, a) e12. That one call replaces the
// clamp-then-acos-then-divide sequence Cl(3,0) needs:
//
//   - It agrees with Cl(3,0)'s formula everywhere. That one computes
//     acos(a) * b/|b|, which for a = cos(phi), b = sin(phi) is
//     |phi| * sign(sin(phi)) = phi on (-pi, pi] -- exactly atan2(b, a).
//   - It is better conditioned near the identity, where acos(a) loses half
//     its precision as a -> 1 and atan2 does not.
//   - It needs no clamp, because it has no domain to leave.
//   - It depends only on the ratio b/a, so a rotor that has drifted off unit
//     length still gives the right angle. (The log of a non-unit complex
//     number also has a real part, ln|r|. A rotor's log is a bivector, so that
//     part is dropped -- as in Cl(3,0).)
//
// The result lies in (-pi, pi], so log(exp(b)) == b for b.xy in that range.
// Plain atan2 would let -pi through when the bivector part is -0.0 and the
// scalar part negative; adding +0.0 turns that -0.0 into +0.0 and changes
// nothing else, so the range really is half-open. geometry.hpp has the detail.

namespace CliffordCore::Cl2
{
    template<typename T>
    /**
     * @brief Computes the logarithm of a rotor in 2D space, resulting in a bivector.
     * @param r The rotor for which to compute the logarithm.
     * @return The resulting bivector, atan2(r.bivector.xy, r.scalar.value) e12,
     *         with its component in (-pi, pi]. The zero rotor gives zero.
     */
    constexpr Bivector<T> log(const Rotor<T>& r) {
        return Bivector<T>(std::atan2(r.bivector.xy + T(0), r.scalar.value));
    }
} // namespace CliffordCore::Cl2
