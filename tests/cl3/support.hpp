// Shared support for the Cl(3,0) test suite: every Cl3 header, the
// helpers typed to Cl3's double types that more than one category file
// uses, and the run_<category>() entry points that main.cpp calls.
//
// Every header is included explicitly rather than through the umbrella,
// because that is what catches a header that fails to pull in its own
// dependencies. Only Cl3 is included: two algebras never meet in one
// translation unit or one binary.

#pragma once

#include "../common/harness.hpp"

#include <type_traits>

#include <cliffordcore/cl3/scalar.hpp>
#include <cliffordcore/cl3/vector.hpp>
#include <cliffordcore/cl3/bivector.hpp>
#include <cliffordcore/cl3/trivector.hpp>
#include <cliffordcore/cl3/multivector.hpp>
#include <cliffordcore/cl3/rotor.hpp>
#include <cliffordcore/cl3/operations/dot_product.hpp>
#include <cliffordcore/cl3/operations/wedge_product.hpp>
#include <cliffordcore/cl3/operations/geometric_product.hpp>
#include <cliffordcore/cl3/operations/norm.hpp>
#include <cliffordcore/cl3/operations/inverse.hpp>
#include <cliffordcore/cl3/operations/reverse.hpp>
#include <cliffordcore/cl3/operations/dual.hpp>
#include <cliffordcore/cl3/operations/exp.hpp>
#include <cliffordcore/cl3/operations/log.hpp>
#include <cliffordcore/cl3/operations/normalize.hpp>
#include <cliffordcore/cl3/operations/sandwich.hpp>
#include <cliffordcore/cl3/operations/addition.hpp>
#include <cliffordcore/cl3/operations/subtraction.hpp>
#include <cliffordcore/cl3/operations/mixed_products.hpp>
#include <cliffordcore/cl3/operations/rotor_construction.hpp>
#include <cliffordcore/cl3/operations/grade.hpp>
#include <cliffordcore/cl3/operations/involutions.hpp>
#include <cliffordcore/cl3/operations/contraction.hpp>
#include <cliffordcore/cl3/operations/geometry.hpp>
#include <cliffordcore/cl3/operations/comparison.hpp>
#include <cliffordcore/cl3/operations/stream.hpp>

namespace cctest {

namespace ga = CliffordCore::Cl3;

using ga::Bivector;
using ga::Multivector;
using ga::Rotor;
using ga::Scalar;
using ga::Trivector;
using ga::Vector;

// Whole-object helpers keep per-component noise out of the test bodies.
inline void check_scalar(const Scalar<double>& s, double value, const std::string& what,
                  double tolerance = 1e-12)
{
    check_close(s.value, value, what + ".value", tolerance);
}

inline void check_vector(const Vector<double>& v, double x, double y, double z, const std::string& what,
                  double tolerance = 1e-12)
{
    check_close(v.x, x, what + ".x", tolerance);
    check_close(v.y, y, what + ".y", tolerance);
    check_close(v.z, z, what + ".z", tolerance);
}

inline void check_bivector(const Bivector<double>& b, double xy, double xz, double yz, const std::string& what,
                    double tolerance = 1e-12)
{
    check_close(b.xy, xy, what + ".xy", tolerance);
    check_close(b.xz, xz, what + ".xz", tolerance);
    check_close(b.yz, yz, what + ".yz", tolerance);
}

inline void check_trivector(const Trivector<double>& t, double e123, const std::string& what,
                     double tolerance = 1e-12)
{
    check_close(t.e123, e123, what + ".e123", tolerance);
}

// Build a multivector from raw components, in grade order.
inline Multivector<double> make_mv(double s, double x, double y, double z,
                             double xy, double xz, double yz, double t)
{
    return Multivector<double>(Scalar<double>(s), Vector<double>(x, y, z),
                                Bivector<double>(xy, xz, yz), Trivector<double>(t));
}

inline void check_mv(const Multivector<double>& m, double s, double x, double y, double z,
              double xy, double xz, double yz, double t, const std::string& what)
{
    check_close(m.scalar.value, s, what + ".scalar");
    check_vector(m.vector, x, y, z, what);
    check_bivector(m.bivector, xy, xz, yz, what);
    check_close(m.trivector.e123, t, what + ".trivector");
}

// Total absolute component difference; 0 means identical.
inline double mv_difference(const Multivector<double>& a, const Multivector<double>& b)
{
    return std::fabs(a.scalar.value - b.scalar.value)
         + std::fabs(a.vector.x - b.vector.x)
         + std::fabs(a.vector.y - b.vector.y)
         + std::fabs(a.vector.z - b.vector.z)
         + std::fabs(a.bivector.xy - b.bivector.xy)
         + std::fabs(a.bivector.xz - b.bivector.xz)
         + std::fabs(a.bivector.yz - b.bivector.yz)
         + std::fabs(a.trivector.e123 - b.trivector.e123);
}

inline Multivector<double> rotor_to_mv(const Rotor<double>& r)
{
    return Multivector<double>(r.scalar, Vector<double>(), r.bivector, Trivector<double>());
}

// One entry point per category file, in the order main.cpp runs them.
void run_types();
void run_arithmetic();
void run_products();
void run_cayley_table();
void run_grades_and_contractions();
void run_unary_ops();
void run_rotations();
void run_geometry();
void run_comparison_and_printing();
void run_properties();
void run_instantiation_sweep();

} // namespace cctest
