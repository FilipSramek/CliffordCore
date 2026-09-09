#pragma once

/**
 * @file format.hpp
 * @brief Shared number formatting for the types' to_string() methods.
 */

#include <limits>
#include <sstream>
#include <string>

namespace CliffordCore
{
namespace detail
{

template<typename T>
/**
 * @brief Formats one component at round-trip precision.
 * @param value The component to format.
 * @return The value as a string that parses back to the same value.
 *
 * Deliberately NOT std::to_string, which is sprintf("%f"): always six decimals
 * and never scientific. That prints 1e-17 as "0.000000" and 1e20 as a
 * 27-character string, and makes float, double and long double look identical.
 * This library routinely produces residuals around 1e-16, so that formatting
 * would be useless for exactly the debugging to_string() exists to support.
 *
 * max_digits10 is the shortest precision that round-trips a floating-point
 * value. It is 0 for integral types, where stream precision is ignored anyway.
 */
inline std::string format_component(T value)
{
    std::ostringstream stream;
    stream.precision(std::numeric_limits<T>::max_digits10);
    stream << value;
    return stream.str();
}

} // namespace detail
} // namespace CliffordCore
