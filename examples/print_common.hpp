#pragma once

// Algebra-independent printing scaffolding for the CliffordCore examples.
//
// Everything here works on plain numbers and strings, so it is shared by the
// Cl(3,0) helpers in print_helpers.hpp and the Cl(3,0,1) ones in
// print_helpers_pga.hpp. Neither algebra's headers are included from here --
// that is what lets a PGA example use the same layout without ever pulling
// Cl3 into the translation unit, which the library design deliberately forbids.
//
// Two choices keep the output reviewable:
//   * every component goes through clean(), which shows values within 1e-12 of
//     zero as plain 0 -- rotations otherwise litter the output with -0.0000 and
//     2.2e-16, which distracts from the point being made;
//   * fixed precision 4 everywhere, formatted into a string, so std::cout's
//     own state is never modified.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace ex {

/// Show near-zero values as exact zero, so rotation output stays readable.
inline double clean(double v)
{
    return std::fabs(v) < 1e-12 ? 0.0 : v;
}

/// Format one component at fixed precision without touching std::cout's state.
inline std::string fmt(double v)
{
    std::ostringstream os;
    os << std::fixed << std::setprecision(4) << clean(v);
    return os.str();
}

/// Format a value in scientific notation, for quantities too small for fmt().
inline std::string sci(double v)
{
    std::ostringstream os;
    os << std::scientific << std::setprecision(3) << v;
    return os.str();
}

inline void heading(const std::string& title)
{
    std::cout << "\n== " << title << " ==\n";
}

/// An indented line of prose explaining why a result looks the way it does.
inline void note(const std::string& text)
{
    std::cout << "   . " << text << "\n";
}

inline void line(const std::string& label, const std::string& value)
{
    std::cout << "  " << std::left << std::setw(30) << label << std::right << value << "\n";
}

} // namespace ex
