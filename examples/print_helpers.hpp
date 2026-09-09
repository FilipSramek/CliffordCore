#pragma once

// Printing scaffolding for the CliffordCore examples.
//
// The library deliberately has no operator<< (see docs/conventions.md, "Known
// gaps"). These helpers exist ONLY so the examples can show their work. They
// live in examples/ rather than include/ so there is no doubt they are not part
// of the API, and they touch nothing but the public component names:
// .value, .x/.y/.z, .xy/.xz/.yz, .e123.
//
// Two deliberate choices keep the output reviewable:
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

#include "cliffordcore.hpp"

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

inline void print(const std::string& label, const CliffordCore::Scalar<double>& s)
{
    line(label, fmt(s.value));
}

inline void print(const std::string& label, const CliffordCore::Vector3<double>& v)
{
    line(label, "(" + fmt(v.x) + ", " + fmt(v.y) + ", " + fmt(v.z) + ")");
}

inline void print(const std::string& label, const CliffordCore::Bivector3<double>& b)
{
    line(label, fmt(b.xy) + "*e12  " + fmt(b.xz) + "*e13  " + fmt(b.yz) + "*e23");
}

inline void print(const std::string& label, const CliffordCore::Trivector3<double>& t)
{
    line(label, fmt(t.e123) + "*e123");
}

inline void print(const std::string& label, const CliffordCore::Rotor3<double>& r)
{
    line(label, fmt(r.scalar.value) + " + " + fmt(r.bivector.xy) + "*e12  "
                + fmt(r.bivector.xz) + "*e13  " + fmt(r.bivector.yz) + "*e23");
}

inline void print(const std::string& label, const CliffordCore::Multivector3<double>& m)
{
    line(label, fmt(m.scalar.value)
                + " | (" + fmt(m.vector.x) + ", " + fmt(m.vector.y) + ", " + fmt(m.vector.z) + ")"
                + " | " + fmt(m.bivector.xy) + " " + fmt(m.bivector.xz) + " " + fmt(m.bivector.yz)
                + " | " + fmt(m.trivector.e123));
}

/// Column header matching the Multivector3 layout above.
inline void multivector_legend()
{
    line("(multivector layout)", "scalar | vector (x,y,z) | bivector xy xz yz | trivector");
}

} // namespace ex
