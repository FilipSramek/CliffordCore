// Shared support for the Cl(2,0) test suite: every Cl2 header, the
// helpers typed to Cl2's double types that more than one category file
// uses, and the run_<category>() entry points that main.cpp calls.
//
// Every header is included explicitly rather than through the umbrella,
// because that is what catches a header that fails to pull in its own
// dependencies. Only Cl2 is included: two algebras never meet in one
// translation unit or one binary.

#pragma once

#include "../common/harness.hpp"

#include <type_traits>
#include <utility>

#include <cliffordcore/cl2/scalar.hpp>
#include <cliffordcore/cl2/vector.hpp>
#include <cliffordcore/cl2/bivector.hpp>
#include <cliffordcore/cl2/multivector.hpp>
#include <cliffordcore/cl2/rotor.hpp>
#include <cliffordcore/cl2/operations/dot_product.hpp>
#include <cliffordcore/cl2/operations/wedge_product.hpp>
#include <cliffordcore/cl2/operations/geometric_product.hpp>
#include <cliffordcore/cl2/operations/mixed_products.hpp>
#include <cliffordcore/cl2/operations/addition.hpp>
#include <cliffordcore/cl2/operations/subtraction.hpp>
#include <cliffordcore/cl2/operations/grade.hpp>
#include <cliffordcore/cl2/operations/contraction.hpp>
#include <cliffordcore/cl2/operations/norm.hpp>
#include <cliffordcore/cl2/operations/normalize.hpp>
#include <cliffordcore/cl2/operations/reverse.hpp>
#include <cliffordcore/cl2/operations/involutions.hpp>
#include <cliffordcore/cl2/operations/inverse.hpp>
#include <cliffordcore/cl2/operations/dual.hpp>
#include <cliffordcore/cl2/operations/comparison.hpp>
#include <cliffordcore/cl2/operations/stream.hpp>
#include <cliffordcore/cl2/operations/exp.hpp>
#include <cliffordcore/cl2/operations/log.hpp>
#include <cliffordcore/cl2/operations/sandwich.hpp>
#include <cliffordcore/cl2/operations/rotor_construction.hpp>
#include <cliffordcore/cl2/operations/geometry.hpp>

namespace cctest {

namespace ga = CliffordCore::Cl2;

using ga::Bivector;
using ga::Multivector;
using ga::Rotor;
using ga::Scalar;
using ga::Vector;

// Whole-object helpers keep per-component noise out of the test bodies.
inline void check_scalar(const Scalar<double>& s, double value, const std::string& what,
                  double tolerance = 1e-12)
{
    check_close(s.value, value, what + ".value", tolerance);
}

inline void check_vector(const Vector<double>& v, double x, double y, const std::string& what,
                  double tolerance = 1e-12)
{
    check_close(v.x, x, what + ".x", tolerance);
    check_close(v.y, y, what + ".y", tolerance);
}

inline void check_bivector(const Bivector<double>& b, double xy, const std::string& what,
                    double tolerance = 1e-12)
{
    check_close(b.xy, xy, what + ".xy", tolerance);
}

inline void check_rotor(const Rotor<double>& r, double s, double xy, const std::string& what,
                 double tolerance = 1e-12)
{
    check_close(r.scalar.value, s, what + ".scalar", tolerance);
    check_close(r.bivector.xy, xy, what + ".bivector.xy", tolerance);
}

// Build a multivector from raw components, in grade order -- the same order as
// the table in docs/cl2.md.
inline Multivector<double> make_mv(double s, double x, double y, double xy)
{
    return Multivector<double>(Scalar<double>(s), Vector<double>(x, y), Bivector<double>(xy));
}

inline void check_mv(const Multivector<double>& m, double s, double x, double y, double xy,
              const std::string& what, double tolerance = 1e-12)
{
    check_close(m.scalar.value, s, what + ".scalar", tolerance);
    check_close(m.vector.x, x, what + ".vector.x", tolerance);
    check_close(m.vector.y, y, what + ".vector.y", tolerance);
    check_close(m.bivector.xy, xy, what + ".bivector.xy", tolerance);
}

inline void check_mv_close(const Multivector<double>& a, const Multivector<double>& b,
                    const std::string& what, double tolerance = 1e-12)
{
    check_mv(a, b.scalar.value, b.vector.x, b.vector.y, b.bivector.xy, what, tolerance);
}

inline void check_vector_close(const Vector<double>& a, const Vector<double>& b,
                        const std::string& what, double tolerance = 1e-12)
{
    check_vector(a, b.x, b.y, what, tolerance);
}

inline void check_rotor_close(const Rotor<double>& a, const Rotor<double>& b,
                       const std::string& what, double tolerance = 1e-12)
{
    check_rotor(a, b.scalar.value, b.bivector.xy, what, tolerance);
}

// Widen a vector into a multivector without reaching into ga::detail, so the
// tests' reference computations use only the public general product.
inline Multivector<double> as_mv(const Vector<double>& v)
{
    return make_mv(0.0, v.x, v.y, 0.0);
}

// The four basis blades, in the library's order, as index lists.
struct BladeIndices
{
    int count;
    int index[2];
};

inline const BladeIndices kBlades[4] = {
    {0, {0, 0}},   // 1
    {1, {1, 0}},   // e1
    {1, {2, 0}},   // e2
    {2, {1, 2}},   // e12
};

inline const char* const kBladeNames[4] = {"1", "e1", "e2", "e12"};

inline Multivector<double> basis_blade(int k)
{
    double c[4] = {0.0, 0.0, 0.0, 0.0};
    c[k] = 1.0;
    return make_mv(c[0], c[1], c[2], c[3]);
}

// Independent reference for the product of two basis blades: concatenate the
// index lists, bubble-sort them counting swaps, then cancel equal neighbours --
// e_i e_i = +1 for both basis vectors, so a cancellation never changes the sign.
// Returns the sign and writes the resulting blade's position.
inline int reference_blade_product(int i, int j, int& result_blade)
{
    int list[4];
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

    int out[2];
    int m = 0;
    for (int k = 0; k < n;) {
        if (k + 1 < n && list[k] == list[k + 1]) {
            k += 2;   // e_i^2 = +1
        } else {
            out[m++] = list[k];
            k += 1;
        }
    }

    result_blade = -1;
    for (int b = 0; b < 4; ++b) {
        if (kBlades[b].count != m) continue;
        bool same = true;
        for (int k = 0; k < m; ++k) same = same && kBlades[b].index[k] == out[k];
        if (same) { result_blade = b; break; }
    }
    return sign;
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
void run_complex_isomorphism();
void run_comparison_and_printing();
void run_properties();
void run_instantiation_sweep();

} // namespace cctest
