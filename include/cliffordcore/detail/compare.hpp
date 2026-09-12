#pragma once

/**
 * @file cliffordcore/detail/compare.hpp
 * @brief Shared tolerance helpers for the approximate comparisons.
 *
 * These are algebra-independent -- they compare raw components, and know
 * nothing about grades or a metric -- so they live beside the shared formatting
 * helper rather than inside any one algebra's namespace. Every algebra's
 * approx_equal builds on them, which keeps the default tolerance consistent
 * across Cl(3,0), PGA and Cl(2,0).
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <cmath>
#include <limits>

namespace CliffordCore
{
namespace detail
{

template<typename T>
/**
 * @brief The default tolerance for approx_equal, scaled to the type.
 * @return One hundred epsilons, which absorbs ordinary rounding without
 *         hiding a real discrepancy.
 */
constexpr T default_tolerance() {
    return std::numeric_limits<T>::epsilon() * T(100);
}

template<typename T>
/**
 * @brief Compares two raw components within a tolerance.
 * @param a The left value.
 * @param b The right value.
 * @param tolerance The largest difference still considered equal.
 * @return True if they differ by no more than the tolerance.
 */
inline bool close(T a, T b, T tolerance) {
    return std::fabs(a - b) <= tolerance;
}

} // namespace detail
} // namespace CliffordCore
