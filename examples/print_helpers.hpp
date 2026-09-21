#pragma once

// Printing scaffolding for the Cl(3,0) CliffordCore examples.
//
// The library now has operator<< (cl3/operations/stream.hpp), which prints at
// round-trip precision. These helpers exist so the examples can show their
// work at a fixed width with near-zero noise suppressed instead. They
// live in examples/ rather than include/ so there is no doubt they are not part
// of the API, and they touch nothing but the public component names:
// .value, .x/.y/.z, .xy/.xz/.yz, .e123.
//
// The formatting itself -- clean(), fmt(), heading(), line() -- is in
// print_common.hpp, shared with the PGA examples. Only the per-type print()
// overloads are here, because those are what differ between algebras.

#include <string>

#include "cliffordcore/cl3.hpp"
#include "print_common.hpp"

namespace ex {

namespace ga = CliffordCore::Cl3;

inline void print(const std::string& label, const ga::Scalar<double>& s)
{
    line(label, fmt(s.value));
}

inline void print(const std::string& label, const ga::Vector<double>& v)
{
    line(label, "(" + fmt(v.x) + ", " + fmt(v.y) + ", " + fmt(v.z) + ")");
}

inline void print(const std::string& label, const ga::Bivector<double>& b)
{
    line(label, fmt(b.xy) + "*e12  " + fmt(b.xz) + "*e13  " + fmt(b.yz) + "*e23");
}

inline void print(const std::string& label, const ga::Trivector<double>& t)
{
    line(label, fmt(t.e123) + "*e123");
}

inline void print(const std::string& label, const ga::Rotor<double>& r)
{
    line(label, fmt(r.scalar.value) + " + " + fmt(r.bivector.xy) + "*e12  "
                + fmt(r.bivector.xz) + "*e13  " + fmt(r.bivector.yz) + "*e23");
}

inline void print(const std::string& label, const ga::Multivector<double>& m)
{
    line(label, fmt(m.scalar.value)
                + " | (" + fmt(m.vector.x) + ", " + fmt(m.vector.y) + ", " + fmt(m.vector.z) + ")"
                + " | " + fmt(m.bivector.xy) + " " + fmt(m.bivector.xz) + " " + fmt(m.bivector.yz)
                + " | " + fmt(m.trivector.e123));
}

/// Column header matching the Multivector layout above.
inline void multivector_legend()
{
    line("(multivector layout)", "scalar | vector (x,y,z) | bivector xy xz yz | trivector");
}

} // namespace ex
