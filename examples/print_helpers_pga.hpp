#pragma once

// Printing scaffolding for the Cl(3,0,1) CliffordCore examples.
//
// The PGA counterpart of print_helpers.hpp. It includes ONLY the PGA umbrella,
// never cl3.hpp: mixing two algebras in one translation unit is the thing the
// namespace-per-algebra design exists to prevent, and an examples helper is no
// exception. The shared formatting lives in print_common.hpp.
//
// Raw components are hard to read as geometry, so the point and line overloads
// show the position and direction alongside them -- which is the whole reason
// primitives.hpp has position() and direction() in the first place.

#include <string>

#include "cliffordcore/pga.hpp"
#include "print_common.hpp"

namespace ex {

namespace pga = CliffordCore::PGA;

inline void print(const std::string& label, const pga::Scalar<double>& s)
{
    line(label, fmt(s.value));
}

/// One term of a plane equation, with its own sign folded into the joining " + "
/// or " - " so a negative coefficient does not read as "+ -3.0000".
inline std::string term(double value, const std::string& symbol)
{
    return (clean(value) < 0.0 ? " - " : " + ") + fmt(std::fabs(clean(value))) + symbol;
}

/// A plane, shown as its equation a x + b y + c z + d = 0.
inline void print(const std::string& label, const pga::Vector<double>& v)
{
    line(label, fmt(v.e1) + "x" + term(v.e2, "y") + term(v.e3, "z") + term(v.e0, "") + " = 0");
}

/// A line, shown as direction plus moment.
inline void print(const std::string& label, const pga::Bivector<double>& b)
{
    const auto d = pga::direction(b);
    line(label, "dir (" + fmt(d[0]) + ", " + fmt(d[1]) + ", " + fmt(d[2]) + ")"
                + "  moment " + fmt(b.e01) + " " + fmt(b.e02) + " " + fmt(b.e03));
}

/// A point, shown as its position; an ideal point as the direction it stands for.
inline void print(const std::string& label, const pga::Trivector<double>& t)
{
    if (pga::is_ideal(t)) {
        const auto d = pga::direction(t);
        line(label, "ideal (" + fmt(d[0]) + ", " + fmt(d[1]) + ", " + fmt(d[2]) + ")");
        return;
    }
    const auto p = pga::position(t);
    line(label, "(" + fmt(p[0]) + ", " + fmt(p[1]) + ", " + fmt(p[2]) + ")"
                + (std::fabs(t.e123 - 1.0) < 1e-12 ? "" : "  weight " + fmt(t.e123)));
}

inline void print(const std::string& label, const pga::Quadvector<double>& q)
{
    line(label, fmt(q.e0123) + "*e0123");
}

inline void print(const std::string& label, const pga::Rotor<double>& r)
{
    line(label, fmt(r.scalar.value) + " + " + fmt(r.e12) + "*e12  "
                + fmt(r.e13) + "*e13  " + fmt(r.e23) + "*e23");
}

inline void print(const std::string& label, const pga::Translator<double>& t)
{
    line(label, fmt(t.scalar.value) + " + " + fmt(t.e01) + "*e01  "
                + fmt(t.e02) + "*e02  " + fmt(t.e03) + "*e03");
}

inline void print(const std::string& label, const pga::Motor<double>& m)
{
    line(label, fmt(m.scalar.value)
                + " | " + fmt(m.bivector.e01) + " " + fmt(m.bivector.e02) + " " + fmt(m.bivector.e03)
                + " | " + fmt(m.bivector.e12) + " " + fmt(m.bivector.e13) + " " + fmt(m.bivector.e23)
                + " | " + fmt(m.quadvector.e0123));
}

/// Column header matching the Motor layout above.
inline void motor_legend()
{
    line("(motor layout)", "scalar | ideal e01 e02 e03 | euclid e12 e13 e23 | e0123");
}

} // namespace ex
