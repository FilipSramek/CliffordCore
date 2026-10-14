// The assertion harness shared by every CliffordCore test suite.
//
// Algebra-neutral on purpose: it includes no algebra, so sharing it never puts
// two algebras in one translation unit -- the same split as
// examples/print_common.hpp. Each suite's support.hpp adds the helpers typed to
// its own algebra on top.
//
// Every check records a pass/fail and the run continues, so a single build
// reports every problem rather than stopping at the first. The counters are
// inline variables, so every translation unit of a suite adds to one count,
// which main() prints.

#pragma once

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace cctest {

inline int g_checks = 0;
inline int g_failures = 0;
inline const char* g_section = "";

inline void section(const char* name)
{
    g_section = name;
}

inline void fail(const std::string& what, const std::string& detail)
{
    ++g_failures;
    std::cerr << "FAIL [" << g_section << "] " << what;
    if (!detail.empty()) {
        std::cerr << " -- " << detail;
    }
    std::cerr << "\n";
}

// Exact boolean check, for things like "this value is identically zero".
inline void check(bool condition, const std::string& what)
{
    ++g_checks;
    if (!condition) {
        fail(what, "");
    }
}

// Floating-point comparison. Never compare doubles with ==; values that are
// mathematically equal can differ in the last bits after arithmetic.
inline void check_close(double actual, double expected, const std::string& what, double tolerance = 1e-12)
{
    ++g_checks;
    if (std::fabs(actual - expected) > tolerance) {
        fail(what, "expected " + std::to_string(expected) + ", got " + std::to_string(actual));
    }
}

inline constexpr double kPi = 3.14159265358979323846;

template<typename T> const char* type_name();
template<> inline const char* type_name<float>() { return "float"; }
template<> inline const char* type_name<double>() { return "double"; }
template<> inline const char* type_name<long double>() { return "long double"; }

// check_close_t allows this many epsilons of T, relative to the size of the
// expected value. A longer computation chain passes a larger `scale`.
inline constexpr int kEpsilons = 64;

// std::to_string prints six fixed decimals, which hides every difference these
// checks care about.
inline std::string describe(long double value)
{
    std::ostringstream os;
    os << std::setprecision(17) << static_cast<double>(value);
    return os.str();
}

template<typename T>
inline void check_close_t(T actual, T expected, const std::string& what, T scale = T(1))
{
    ++g_checks;
    const T tolerance = T(kEpsilons) * std::numeric_limits<T>::epsilon()
                      * std::max(T(1), std::fabs(expected)) * scale;
    // Negated so that a NaN fails instead of slipping through.
    if (!(std::fabs(actual - expected) <= tolerance)) {
        fail(what, std::string(type_name<T>()) + ": expected " + describe(expected)
                   + ", got " + describe(actual) + ", tolerance " + describe(tolerance));
    }
}

} // namespace cctest
