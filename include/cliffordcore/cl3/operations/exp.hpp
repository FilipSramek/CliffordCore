#pragma once

/**
 * @file cliffordcore/cl3/operations/exp.hpp
 * @brief Exponential of a bivector, producing a rotor.
 *
 * exp(theta * B) is the rotor rotating by 2 * theta in the plane B, so a
 * quarter turn is exp(pi/4 * B). The unit bivector squares to -1, which
 * makes the series collapse to cos + sin exactly as it does for complex
 * numbers.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
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
     * @brief Computes the exponential of a bivector in 3D space, resulting in a rotor.
     * @param b The bivector for which to compute the exponential.
     * @return The resulting rotor from the exponential of the bivector b.
     */
    constexpr Rotor<T> exp(const Bivector<T>& b) {
        Scalar<T> norm_bivector = norm(b);

        if (norm_bivector.value == 0) { // Handle the case when the bivector is zero
            return Rotor<T>(Scalar<T>(1), Bivector<T>(0, 0, 0));
        }

        Scalar<T> cos_norm = Scalar<T>(std::cos(norm_bivector.value));

        Bivector<T> sin_norm = Bivector<T>(
            std::sin(norm_bivector.value) * (b.xy / norm_bivector.value),
            std::sin(norm_bivector.value) * (b.xz / norm_bivector.value),
            std::sin(norm_bivector.value) * (b.yz / norm_bivector.value)
        );

        return Rotor<T>(cos_norm, sin_norm);
    }
} // namespace CliffordCore::Cl3
