// Shared support for the Cl(3,0,1) test suite: every PGA header, the
// helpers typed to PGA's double types that more than one category file
// uses, and the run_<category>() entry points that main.cpp calls.
//
// Every header is included explicitly rather than through the umbrella,
// because that is what catches a header that fails to pull in its own
// dependencies. Only PGA is included: two algebras never meet in one
// translation unit or one binary.

#pragma once

#include "../common/harness.hpp"

#include <array>
#include <type_traits>

#include <cliffordcore/pga/scalar.hpp>
#include <cliffordcore/pga/vector.hpp>
#include <cliffordcore/pga/bivector.hpp>
#include <cliffordcore/pga/trivector.hpp>
#include <cliffordcore/pga/quadvector.hpp>
#include <cliffordcore/pga/multivector.hpp>
#include <cliffordcore/pga/rotor.hpp>
#include <cliffordcore/pga/translator.hpp>
#include <cliffordcore/pga/motor.hpp>
#include <cliffordcore/pga/operations/dot_product.hpp>
#include <cliffordcore/pga/operations/wedge_product.hpp>
#include <cliffordcore/pga/operations/geometric_product.hpp>
#include <cliffordcore/pga/operations/mixed_products.hpp>
#include <cliffordcore/pga/operations/addition.hpp>
#include <cliffordcore/pga/operations/subtraction.hpp>
#include <cliffordcore/pga/operations/grade.hpp>
#include <cliffordcore/pga/operations/norm.hpp>
#include <cliffordcore/pga/operations/normalize.hpp>
#include <cliffordcore/pga/operations/reverse.hpp>
#include <cliffordcore/pga/operations/involutions.hpp>
#include <cliffordcore/pga/operations/inverse.hpp>
#include <cliffordcore/pga/operations/dual.hpp>
#include <cliffordcore/pga/operations/regressive_product.hpp>
#include <cliffordcore/pga/operations/contraction.hpp>
#include <cliffordcore/pga/operations/comparison.hpp>
#include <cliffordcore/pga/operations/stream.hpp>
#include <cliffordcore/pga/operations/exp.hpp>
#include <cliffordcore/pga/operations/log.hpp>
#include <cliffordcore/pga/operations/sandwich.hpp>
#include <cliffordcore/pga/operations/primitives.hpp>
#include <cliffordcore/pga/operations/motor_construction.hpp>
#include <cliffordcore/pga/operations/geometry.hpp>

namespace cctest {

namespace ga = CliffordCore::PGA;

using ga::Bivector;
using ga::Motor;
using ga::Multivector;
using ga::Quadvector;
using ga::Rotor;
using ga::Scalar;
using ga::Translator;
using ga::Trivector;
using ga::Vector;

// Whole-object helpers keep per-component noise out of the test bodies.
inline void check_scalar(const Scalar<double>& s, double value, const std::string& what,
                  double tolerance = 1e-12)
{
    check_close(s.value, value, what + ".value", tolerance);
}

inline void check_vector(const Vector<double>& v, double e0, double e1, double e2, double e3,
                  const std::string& what, double tolerance = 1e-12)
{
    check_close(v.e0, e0, what + ".e0", tolerance);
    check_close(v.e1, e1, what + ".e1", tolerance);
    check_close(v.e2, e2, what + ".e2", tolerance);
    check_close(v.e3, e3, what + ".e3", tolerance);
}

inline void check_bivector(const Bivector<double>& b, double e01, double e02, double e03,
                    double e12, double e13, double e23, const std::string& what,
                    double tolerance = 1e-12)
{
    check_close(b.e01, e01, what + ".e01", tolerance);
    check_close(b.e02, e02, what + ".e02", tolerance);
    check_close(b.e03, e03, what + ".e03", tolerance);
    check_close(b.e12, e12, what + ".e12", tolerance);
    check_close(b.e13, e13, what + ".e13", tolerance);
    check_close(b.e23, e23, what + ".e23", tolerance);
}

inline void check_trivector(const Trivector<double>& t, double e012, double e013, double e023, double e123,
                     const std::string& what, double tolerance = 1e-12)
{
    check_close(t.e012, e012, what + ".e012", tolerance);
    check_close(t.e013, e013, what + ".e013", tolerance);
    check_close(t.e023, e023, what + ".e023", tolerance);
    check_close(t.e123, e123, what + ".e123", tolerance);
}

inline void check_quadvector(const Quadvector<double>& q, double e0123, const std::string& what,
                      double tolerance = 1e-12)
{
    check_close(q.e0123, e0123, what + ".e0123", tolerance);
}

inline void check_rotor(const Rotor<double>& r, double s, double e12, double e13, double e23,
                 const std::string& what, double tolerance = 1e-12)
{
    check_close(r.scalar.value, s, what + ".scalar", tolerance);
    check_close(r.e12, e12, what + ".e12", tolerance);
    check_close(r.e13, e13, what + ".e13", tolerance);
    check_close(r.e23, e23, what + ".e23", tolerance);
}

inline void check_translator(const Translator<double>& t, double s, double e01, double e02, double e03,
                      const std::string& what, double tolerance = 1e-12)
{
    check_close(t.scalar.value, s, what + ".scalar", tolerance);
    check_close(t.e01, e01, what + ".e01", tolerance);
    check_close(t.e02, e02, what + ".e02", tolerance);
    check_close(t.e03, e03, what + ".e03", tolerance);
}

// Build a multivector from raw components, in grade order and lexicographic
// order within each grade -- the same order as the table in docs/pga.md.
inline Multivector<double> make_mv(double s,
                            double e0, double e1, double e2, double e3,
                            double e01, double e02, double e03, double e12, double e13, double e23,
                            double e012, double e013, double e023, double e123,
                            double e0123)
{
    return Multivector<double>(
        Scalar<double>(s),
        Vector<double>(e0, e1, e2, e3),
        Bivector<double>(e01, e02, e03, e12, e13, e23),
        Trivector<double>(e012, e013, e023, e123),
        Quadvector<double>(e0123));
}

// Total absolute component difference; 0 means identical.
inline double mv_difference(const Multivector<double>& a, const Multivector<double>& b)
{
    return std::fabs(a.scalar.value - b.scalar.value)
         + std::fabs(a.vector.e0 - b.vector.e0)
         + std::fabs(a.vector.e1 - b.vector.e1)
         + std::fabs(a.vector.e2 - b.vector.e2)
         + std::fabs(a.vector.e3 - b.vector.e3)
         + std::fabs(a.bivector.e01 - b.bivector.e01)
         + std::fabs(a.bivector.e02 - b.bivector.e02)
         + std::fabs(a.bivector.e03 - b.bivector.e03)
         + std::fabs(a.bivector.e12 - b.bivector.e12)
         + std::fabs(a.bivector.e13 - b.bivector.e13)
         + std::fabs(a.bivector.e23 - b.bivector.e23)
         + std::fabs(a.trivector.e012 - b.trivector.e012)
         + std::fabs(a.trivector.e013 - b.trivector.e013)
         + std::fabs(a.trivector.e023 - b.trivector.e023)
         + std::fabs(a.trivector.e123 - b.trivector.e123)
         + std::fabs(a.quadvector.e0123 - b.quadvector.e0123);
}

// Sixteen components is too many to spell out per call, so whole-multivector
// expectations are written as another multivector and compared as a whole.
inline void check_mv_close(const Multivector<double>& actual, const Multivector<double>& expected,
                    const std::string& what, double tolerance = 1e-12)
{
    ++g_checks;
    const double diff = mv_difference(actual, expected);
    if (diff > tolerance) {
        fail(what, "total component difference " + std::to_string(diff)
                   + "\n    got      " + actual.to_string()
                   + "\n    expected " + expected.to_string());
    }
}

// The sixteen basis blades in lexicographic order, as index lists. Index 0 is
// the scalar (empty list). This order is the one make_mv() takes.
struct BladeIndices {
    int count;
    int index[4];
};

inline const BladeIndices kBlades[16] = {
    {0, {0, 0, 0, 0}},
    {1, {0, 0, 0, 0}}, {1, {1, 0, 0, 0}}, {1, {2, 0, 0, 0}}, {1, {3, 0, 0, 0}},
    {2, {0, 1, 0, 0}}, {2, {0, 2, 0, 0}}, {2, {0, 3, 0, 0}},
    {2, {1, 2, 0, 0}}, {2, {1, 3, 0, 0}}, {2, {2, 3, 0, 0}},
    {3, {0, 1, 2, 0}}, {3, {0, 1, 3, 0}}, {3, {0, 2, 3, 0}}, {3, {1, 2, 3, 0}},
    {4, {0, 1, 2, 3}},
};

inline const char* const kBladeNames[16] = {
    "1", "e0", "e1", "e2", "e3", "e01", "e02", "e03", "e12", "e13", "e23",
    "e012", "e013", "e023", "e123", "e0123",
};

// A multivector with a single basis blade set to one.
inline Multivector<double> basis_blade(int which)
{
    double c[16] = {0};
    c[which] = 1.0;
    return make_mv(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8],
                   c[9], c[10], c[11], c[12], c[13], c[14], c[15]);
}

// Independent reference for the product of two basis blades: concatenate the
// index lists, bubble-sort them counting swaps (each swap flips the sign),
// then cancel equal neighbours -- e0 e0 = 0 kills the term, e_i e_i = +1 for
// the others. Returns the sign (0 if the product vanishes) and the index of
// the resulting blade. This is the definition of the algebra; the 16x16 table
// in geometric_product.hpp must reproduce it exactly.
inline int reference_blade_product(int i, int j, int& result_blade)
{
    int list[8];
    int n = 0;
    for (int k = 0; k < kBlades[i].count; ++k) list[n++] = kBlades[i].index[k];
    for (int k = 0; k < kBlades[j].count; ++k) list[n++] = kBlades[j].index[k];

    int sign = 1;
    bool changed = true;
    while (changed) {
        changed = false;
        for (int k = 0; k + 1 < n; ++k) {
            if (list[k] > list[k + 1]) {
                const int tmp = list[k]; list[k] = list[k + 1]; list[k + 1] = tmp;
                sign = -sign;
                changed = true;
            }
        }
    }

    int out[4];
    int m = 0;
    for (int k = 0; k < n;) {
        if (k + 1 < n && list[k] == list[k + 1]) {
            if (list[k] == 0) sign = 0;   // e0^2 = 0
            k += 2;
        } else {
            out[m++] = list[k];
            k += 1;
        }
    }

    result_blade = -1;
    for (int b = 0; b < 16; ++b) {
        if (kBlades[b].count != m) continue;
        bool same = true;
        for (int k = 0; k < m; ++k) same = same && kBlades[b].index[k] == out[k];
        if (same) { result_blade = b; break; }
    }
    return sign;
}

// The primitives under short local names, so the tests read as geometry.
// These are the real factories from primitives.hpp.
inline Trivector<double> test_point(double x, double y, double z)
{
    return ga::point(x, y, z);
}

inline Vector<double> test_plane(double a, double b, double c, double d)
{
    return ga::plane(a, b, c, d);
}

// One entry point per category file, in the order main.cpp runs them.
void run_types();
void run_arithmetic();
void run_products();
void run_cayley_table();
void run_incidence();
void run_grades_and_contractions();
void run_unary_ops();
void run_motions();
void run_geometry();
void run_comparison_and_printing();
void run_properties();
void run_instantiation_sweep();

} // namespace cctest
