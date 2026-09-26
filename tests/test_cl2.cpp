// CliffordCore test suite for Cl(2,0), namespace CliffordCore::Cl2.
//
// Build and run:
//   g++ -std=c++17 -Iinclude tests/test_cl2.cpp -o build/test_cl2 && ./build/test_cl2.exe
//   ...or just ./build.sh, which builds every tests/*.cpp
//
// Structure mirrors tests/test_core.cpp and tests/test_pga.cpp: a small
// assertion harness, then one section per type and per operation. Every check
// records a pass/fail and the run continues, so a single build reports every
// problem rather than stopping at the first. This file deliberately includes
// each header explicitly rather than the umbrella, because that is what catches
// a header that fails to pull in its own dependencies.
//
// The sign conventions asserted here are the ones documented in docs/cl2.md;
// every one of them was verified against an independent model of the algebra
// before being written down, so a failure means the code, not the expectation.
//
// This file -- and only this file -- includes <complex>. The even subalgebra of
// Cl(2,0) is the complex numbers, so std::complex is an independent reference
// implementation of every rotor operation, and test_complex_isomorphism checks
// the library against it. No header may take that dependency.

#include <cmath>
#include <complex>
#include <iostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>

#include "../include/cliffordcore/cl2/scalar.hpp"
#include "../include/cliffordcore/cl2/vector.hpp"
#include "../include/cliffordcore/cl2/bivector.hpp"
#include "../include/cliffordcore/cl2/multivector.hpp"
#include "../include/cliffordcore/cl2/rotor.hpp"
#include "../include/cliffordcore/cl2/operations/dot_product.hpp"
#include "../include/cliffordcore/cl2/operations/wedge_product.hpp"
#include "../include/cliffordcore/cl2/operations/geometric_product.hpp"
#include "../include/cliffordcore/cl2/operations/mixed_products.hpp"
#include "../include/cliffordcore/cl2/operations/addition.hpp"
#include "../include/cliffordcore/cl2/operations/subtraction.hpp"
#include "../include/cliffordcore/cl2/operations/grade.hpp"
#include "../include/cliffordcore/cl2/operations/contraction.hpp"
#include "../include/cliffordcore/cl2/operations/norm.hpp"
#include "../include/cliffordcore/cl2/operations/normalize.hpp"
#include "../include/cliffordcore/cl2/operations/reverse.hpp"
#include "../include/cliffordcore/cl2/operations/involutions.hpp"
#include "../include/cliffordcore/cl2/operations/inverse.hpp"
#include "../include/cliffordcore/cl2/operations/dual.hpp"
#include "../include/cliffordcore/cl2/operations/comparison.hpp"
#include "../include/cliffordcore/cl2/operations/stream.hpp"
#include "../include/cliffordcore/cl2/operations/exp.hpp"
#include "../include/cliffordcore/cl2/operations/log.hpp"
#include "../include/cliffordcore/cl2/operations/sandwich.hpp"
#include "../include/cliffordcore/cl2/operations/rotor_construction.hpp"
#include "../include/cliffordcore/cl2/operations/geometry.hpp"

namespace ga = CliffordCore::Cl2;

using ga::Bivector;
using ga::Multivector;
using ga::Rotor;
using ga::Scalar;
using ga::Vector;

namespace {

// ---------------------------------------------------------------------------
// Assertion harness
// ---------------------------------------------------------------------------

int g_checks = 0;
int g_failures = 0;
const char* g_section = "";

void section(const char* name)
{
    g_section = name;
}

void fail(const std::string& what, const std::string& detail)
{
    ++g_failures;
    std::cerr << "FAIL [" << g_section << "] " << what;
    if (!detail.empty()) {
        std::cerr << " -- " << detail;
    }
    std::cerr << "\n";
}

// Exact boolean check, for things like "this value is identically zero".
void check(bool condition, const std::string& what)
{
    ++g_checks;
    if (!condition) {
        fail(what, "");
    }
}

// Floating-point comparison. Never compare doubles with ==; values that are
// mathematically equal can differ in the last bits after arithmetic.
void check_close(double actual, double expected, const std::string& what, double tolerance = 1e-12)
{
    ++g_checks;
    if (std::fabs(actual - expected) > tolerance) {
        fail(what, "expected " + std::to_string(expected) + ", got " + std::to_string(actual));
    }
}

// Whole-object helpers keep per-component noise out of the test bodies.
void check_scalar(const Scalar<double>& s, double value, const std::string& what,
                  double tolerance = 1e-12)
{
    check_close(s.value, value, what + ".value", tolerance);
}

void check_vector(const Vector<double>& v, double x, double y, const std::string& what,
                  double tolerance = 1e-12)
{
    check_close(v.x, x, what + ".x", tolerance);
    check_close(v.y, y, what + ".y", tolerance);
}

void check_bivector(const Bivector<double>& b, double xy, const std::string& what,
                    double tolerance = 1e-12)
{
    check_close(b.xy, xy, what + ".xy", tolerance);
}

void check_rotor(const Rotor<double>& r, double s, double xy, const std::string& what,
                 double tolerance = 1e-12)
{
    check_close(r.scalar.value, s, what + ".scalar", tolerance);
    check_close(r.bivector.xy, xy, what + ".bivector.xy", tolerance);
}

// Build a multivector from raw components, in grade order -- the same order as
// the table in docs/cl2.md.
Multivector<double> make_mv(double s, double x, double y, double xy)
{
    return Multivector<double>(Scalar<double>(s), Vector<double>(x, y), Bivector<double>(xy));
}

// Total absolute component difference; 0 means identical.
double mv_difference(const Multivector<double>& a, const Multivector<double>& b)
{
    return std::fabs(a.scalar.value - b.scalar.value)
         + std::fabs(a.vector.x - b.vector.x)
         + std::fabs(a.vector.y - b.vector.y)
         + std::fabs(a.bivector.xy - b.bivector.xy);
}

void check_mv(const Multivector<double>& m, double s, double x, double y, double xy,
              const std::string& what, double tolerance = 1e-12)
{
    check_close(m.scalar.value, s, what + ".scalar", tolerance);
    check_close(m.vector.x, x, what + ".vector.x", tolerance);
    check_close(m.vector.y, y, what + ".vector.y", tolerance);
    check_close(m.bivector.xy, xy, what + ".bivector.xy", tolerance);
}

void check_mv_close(const Multivector<double>& a, const Multivector<double>& b,
                    const std::string& what, double tolerance = 1e-12)
{
    check_mv(a, b.scalar.value, b.vector.x, b.vector.y, b.bivector.xy, what, tolerance);
}

void check_vector_close(const Vector<double>& a, const Vector<double>& b,
                        const std::string& what, double tolerance = 1e-12)
{
    check_vector(a, b.x, b.y, what, tolerance);
}

void check_rotor_close(const Rotor<double>& a, const Rotor<double>& b,
                       const std::string& what, double tolerance = 1e-12)
{
    check_rotor(a, b.scalar.value, b.bivector.xy, what, tolerance);
}

constexpr double kPi = 3.14159265358979323846;

// Widen a vector into a multivector without reaching into ga::detail, so the
// reference computations below use only the public general product.
Multivector<double> as_mv(const Vector<double>& v)
{
    return make_mv(0.0, v.x, v.y, 0.0);
}

// Name an operand pair's product type, for the return-type table.
template<typename A, typename B>
using product_t = decltype(std::declval<A>() * std::declval<B>());

template<typename A, typename B>
using sum_t = decltype(std::declval<A>() + std::declval<B>());

template<typename A, typename B>
using difference_t = decltype(std::declval<A>() - std::declval<B>());

// Detects whether a ^ b compiles, so a deliberately missing overload can be
// asserted missing rather than merely not tested.
template<typename A, typename B, typename = void>
struct has_wedge : std::false_type {};

template<typename A, typename B>
struct has_wedge<A, B, decltype(void(std::declval<A>() ^ std::declval<B>()))> : std::true_type {};

// ---------------------------------------------------------------------------
// Construction and defaults
// ---------------------------------------------------------------------------

void test_construction()
{
    section("construction");

    check_scalar(Scalar<double>(), 0.0, "default Scalar");
    check_vector(Vector<double>(), 0.0, 0.0, "default Vector");
    check_bivector(Bivector<double>(), 0.0, "default Bivector");
    check_rotor(Rotor<double>(), 0.0, 0.0, "default Rotor");
    check_mv(Multivector<double>(), 0.0, 0.0, 0.0, 0.0, "default Multivector");

    check_scalar(Scalar<double>(2.5), 2.5, "Scalar(2.5)");
    check_vector(Vector<double>(3.0, -4.0), 3.0, -4.0, "Vector(3, -4)");
    check_bivector(Bivector<double>(7.0), 7.0, "Bivector(7)");
    check_rotor(Rotor<double>(Scalar<double>(1.0), Bivector<double>(2.0)), 1.0, 2.0, "Rotor(1, 2)");
    check_mv(make_mv(1.0, 2.0, 3.0, 4.0), 1.0, 2.0, 3.0, 4.0, "make_mv");

    // magnitude() on the blade types.
    check_close(Vector<double>(3.0, 4.0).magnitude().value, 5.0, "Vector(3,4).magnitude");
    check_close(Bivector<double>(-6.0).magnitude().value, 6.0, "Bivector(-6).magnitude");
}

// ---------------------------------------------------------------------------
// Same-type arithmetic
// ---------------------------------------------------------------------------

void test_same_type_arithmetic()
{
    section("same-type arithmetic");

    const Scalar<double> s1(3.0), s2(4.0);
    check_scalar(s1 + s2, 7.0, "Scalar +");
    check_scalar(s1 - s2, -1.0, "Scalar -");
    check_scalar(s1 * s2, 12.0, "Scalar *");
    check_scalar(s1 / s2, 0.75, "Scalar /");
    check_scalar(-s1, -3.0, "Scalar unary -");

    const Vector<double> v1(1.0, 2.0), v2(4.0, 8.0);
    check_vector(v1 + v2, 5.0, 10.0, "Vector +");
    check_vector(v1 - v2, -3.0, -6.0, "Vector -");
    check_vector(-v1, -1.0, -2.0, "Vector unary -");

    const Bivector<double> b1(3.0), b2(5.0);
    check_bivector(b1 + b2, 8.0, "Bivector +");
    check_bivector(b1 - b2, -2.0, "Bivector -");
    check_bivector(-b1, -3.0, "Bivector unary -");

    const Rotor<double> r1(Scalar<double>(1.0), Bivector<double>(2.0));
    const Rotor<double> r2(Scalar<double>(3.0), Bivector<double>(5.0));
    check_rotor(r1 + r2, 4.0, 7.0, "Rotor +");
    check_rotor(r1 - r2, -2.0, -3.0, "Rotor -");

    const Multivector<double> m1 = make_mv(1.0, 2.0, 3.0, 4.0);
    const Multivector<double> m2 = make_mv(10.0, 20.0, 30.0, 40.0);
    check_mv(m1 + m2, 11.0, 22.0, 33.0, 44.0, "Multivector +");
    check_mv(m2 - m1, 9.0, 18.0, 27.0, 36.0, "Multivector -");
    check_mv(-m1, -1.0, -2.0, -3.0, -4.0, "Multivector unary -");
}

// ---------------------------------------------------------------------------
// The pseudoscalar squares to -1
// ---------------------------------------------------------------------------

void test_bivector_is_the_pseudoscalar()
{
    section("bivector as pseudoscalar");

    // e12 * e12 = -1, exactly as Cl(3,0)'s e123 does, so two bivectors multiply
    // to a Scalar rather than to a Bivector. This is the type-level fact that
    // separates Cl(2,0) from Cl(3,0), where Bivector * Bivector is mixed grade.
    const Bivector<double> I(1.0);
    check_scalar(I * I, -1.0, "e12 * e12");

    const Bivector<double> a(3.0), b(2.0);
    check_scalar(a * b, -6.0, "(3 e12)(2 e12)");
    static_assert(std::is_same<decltype(a * b), Scalar<double>>::value,
                  "Bivector * Bivector must return Scalar");

    // Division cancels both minus signs, so it is plain a/b.
    check_scalar(a / b, 1.5, "(3 e12)/(2 e12)");
    static_assert(std::is_same<decltype(a / b), Scalar<double>>::value,
                  "Bivector / Bivector must return Scalar");

    // Scaling a bivector still yields a Bivector, so b * 2.0 and b * b2
    // deliberately differ in return type. These are the exact-match raw-T
    // overloads that keep that unambiguous -- without them the line below is a
    // compile error, because a raw double can reach either operator* by one
    // user-defined conversion.
    check_bivector(a * 2.0, 6.0, "Bivector * raw T");
    check_bivector(2.0 * a, 6.0, "raw T * Bivector");
    check_bivector(a / 2.0, 1.5, "Bivector / raw T");
    static_assert(std::is_same<decltype(a * 2.0), Bivector<double>>::value,
                  "Bivector * T must return Bivector");
}

// ---------------------------------------------------------------------------
// Scalar mixing, on either side
// ---------------------------------------------------------------------------

void test_scalar_products()
{
    section("scalar products");

    const Scalar<double> two(2.0);

    const Vector<double> v(1.0, -3.0);
    check_vector(two * v, 2.0, -6.0, "Scalar * Vector");
    check_vector(v * two, 2.0, -6.0, "Vector * Scalar");
    check_vector(v / two, 0.5, -1.5, "Vector / Scalar");
    check_vector(2.0 * v, 2.0, -6.0, "raw T * Vector");

    const Bivector<double> b(5.0);
    check_bivector(two * b, 10.0, "Scalar * Bivector");
    check_bivector(b * two, 10.0, "Bivector * Scalar");
    check_bivector(b / two, 2.5, "Bivector / Scalar");

    const Rotor<double> r(Scalar<double>(1.0), Bivector<double>(4.0));
    check_rotor(two * r, 2.0, 8.0, "Scalar * Rotor");
    check_rotor(r * two, 2.0, 8.0, "Rotor * Scalar");
    check_rotor(r / two, 0.5, 2.0, "Rotor / Scalar");
    check_rotor(2.0 * r, 2.0, 8.0, "raw T * Rotor");

    const Multivector<double> m = make_mv(1.0, 2.0, 3.0, 4.0);
    check_mv(two * m, 2.0, 4.0, 6.0, 8.0, "Scalar * Multivector");
    check_mv(m * two, 2.0, 4.0, 6.0, 8.0, "Multivector * Scalar");
    check_mv(m / two, 0.5, 1.0, 1.5, 2.0, "Multivector / Scalar");
    check_mv(2.0 * m, 2.0, 4.0, 6.0, 8.0, "raw T * Multivector");

    // Return types: scaling never changes the grade.
    static_assert(std::is_same<decltype(two * v), Vector<double>>::value, "Scalar * Vector");
    static_assert(std::is_same<decltype(two * b), Bivector<double>>::value, "Scalar * Bivector");
    static_assert(std::is_same<decltype(two * r), Rotor<double>>::value, "Scalar * Rotor");
    static_assert(std::is_same<decltype(two * m), Multivector<double>>::value, "Scalar * Multivector");
}

// ---------------------------------------------------------------------------
// Compound assignment
// ---------------------------------------------------------------------------

void test_compound_assignment()
{
    section("compound assignment");

    Scalar<double> s(10.0);
    s += Scalar<double>(5.0);  check_scalar(s, 15.0, "Scalar +=");
    s -= Scalar<double>(3.0);  check_scalar(s, 12.0, "Scalar -=");
    s *= Scalar<double>(2.0);  check_scalar(s, 24.0, "Scalar *=");
    s /= Scalar<double>(4.0);  check_scalar(s, 6.0, "Scalar /=");

    Vector<double> v(1.0, 2.0);
    v += Vector<double>(3.0, 4.0);  check_vector(v, 4.0, 6.0, "Vector +=");
    v -= Vector<double>(1.0, 1.0);  check_vector(v, 3.0, 5.0, "Vector -=");
    v *= 2.0;                       check_vector(v, 6.0, 10.0, "Vector *= raw T");
    v /= 2.0;                       check_vector(v, 3.0, 5.0, "Vector /= raw T");
    v *= Scalar<double>(3.0);       check_vector(v, 9.0, 15.0, "Vector *= Scalar");
    v /= Scalar<double>(3.0);       check_vector(v, 3.0, 5.0, "Vector /= Scalar");

    Bivector<double> b(4.0);
    b += Bivector<double>(2.0);  check_bivector(b, 6.0, "Bivector +=");
    b -= Bivector<double>(1.0);  check_bivector(b, 5.0, "Bivector -=");
    b *= 2.0;                    check_bivector(b, 10.0, "Bivector *= raw T");
    b /= 5.0;                    check_bivector(b, 2.0, "Bivector /= raw T");
    b *= Scalar<double>(3.0);    check_bivector(b, 6.0, "Bivector *= Scalar");
    b /= Scalar<double>(2.0);    check_bivector(b, 3.0, "Bivector /= Scalar");

    Rotor<double> r(Scalar<double>(1.0), Bivector<double>(2.0));
    r += Rotor<double>(Scalar<double>(3.0), Bivector<double>(4.0));
    check_rotor(r, 4.0, 6.0, "Rotor +=");
    r -= Rotor<double>(Scalar<double>(1.0), Bivector<double>(1.0));
    check_rotor(r, 3.0, 5.0, "Rotor -=");
    r *= Scalar<double>(2.0);  check_rotor(r, 6.0, 10.0, "Rotor *=");
    r /= Scalar<double>(2.0);  check_rotor(r, 3.0, 5.0, "Rotor /=");

    Multivector<double> m = make_mv(1.0, 2.0, 3.0, 4.0);
    m += make_mv(1.0, 1.0, 1.0, 1.0);  check_mv(m, 2.0, 3.0, 4.0, 5.0, "Multivector +=");
    m -= make_mv(1.0, 1.0, 1.0, 1.0);  check_mv(m, 1.0, 2.0, 3.0, 4.0, "Multivector -=");
    m *= Scalar<double>(2.0);          check_mv(m, 2.0, 4.0, 6.0, 8.0, "Multivector *=");
    m /= Scalar<double>(2.0);          check_mv(m, 1.0, 2.0, 3.0, 4.0, "Multivector /=");
}

// ---------------------------------------------------------------------------
// Conversions between the even types and the general one
// ---------------------------------------------------------------------------

void test_conversions()
{
    section("conversions");

    // Widening a rotor is lossless, so it is implicit.
    static_assert(std::is_convertible<Rotor<double>, Multivector<double>>::value,
                  "Rotor must widen implicitly to Multivector");

    // Narrowing discards grade 1, so it must be spelled out. This is the trap
    // Cl(3,0) hit: while the constructor was implicit, `rotor + multivector`
    // silently threw grades away.
    static_assert(!std::is_convertible<Multivector<double>, Rotor<double>>::value,
                  "Multivector must NOT convert implicitly to Rotor");
    static_assert(std::is_constructible<Rotor<double>, Multivector<double>>::value,
                  "Rotor must be explicitly constructible from Multivector");

    const Rotor<double> r(Scalar<double>(2.0), Bivector<double>(3.0));
    const Multivector<double> widened = r;
    check_mv(widened, 2.0, 0.0, 0.0, 3.0, "Rotor widened to Multivector");

    const Multivector<double> m = make_mv(5.0, 6.0, 7.0, 8.0);
    const Rotor<double> narrowed(m);
    check_rotor(narrowed, 5.0, 8.0, "Multivector narrowed to Rotor");

    // The aliases name the same types, not new ones.
    static_assert(std::is_same<ga::Pseudoscalar<double>, Bivector<double>>::value,
                  "Pseudoscalar must alias Bivector");
    static_assert(std::is_same<ga::Complex<double>, Rotor<double>>::value,
                  "Complex must alias Rotor");
}

// ---------------------------------------------------------------------------
// String formatting
// ---------------------------------------------------------------------------

void test_to_string()
{
    section("to_string");

    check(Scalar<double>(1.5).to_string() == "1.5", "Scalar::to_string");
    check(Vector<double>(1.0, 2.0).to_string() == "1*e1 + 2*e2", "Vector::to_string");
    check(Bivector<double>(3.0).to_string() == "3*e12", "Bivector::to_string");
    check(Rotor<double>(Scalar<double>(1.0), Bivector<double>(2.0)).to_string()
              == "(1) + (2*e12)",
          "Rotor::to_string");
    check(make_mv(1.0, 2.0, 3.0, 4.0).to_string() == "(1) + (2*e1 + 3*e2) + (4*e12)",
          "Multivector::to_string");

    // Round-trip precision, not std::to_string's fixed six decimals: a 1e-16
    // residual must not print as 0.000000.
    check(Scalar<double>(1e-16).to_string() != "0.000000", "to_string keeps small values");
}

// ---------------------------------------------------------------------------
// The multiplication table
// ---------------------------------------------------------------------------

// The four basis blades, in the library's order, as index lists.
struct BladeIndices
{
    int count;
    int index[2];
};

const BladeIndices kBlades[4] = {
    {0, {0, 0}},   // 1
    {1, {1, 0}},   // e1
    {1, {2, 0}},   // e2
    {2, {1, 2}},   // e12
};

const char* const kBladeNames[4] = {"1", "e1", "e2", "e12"};

Multivector<double> basis_blade(int k)
{
    double c[4] = {0.0, 0.0, 0.0, 0.0};
    c[k] = 1.0;
    return make_mv(c[0], c[1], c[2], c[3]);
}

// Independent reference for the product of two basis blades: concatenate the
// index lists, bubble-sort them counting swaps, then cancel equal neighbours --
// e_i e_i = +1 for both basis vectors, so a cancellation never changes the sign.
// Returns the sign and writes the resulting blade's position.
int reference_blade_product(int i, int j, int& result_blade)
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

void test_basis_multiplication_table()
{
    section("basis multiplication table");

    const Multivector<double> one = basis_blade(0);
    const Multivector<double> e1 = basis_blade(1);
    const Multivector<double> e2 = basis_blade(2);
    const Multivector<double> e12 = basis_blade(3);

    check_mv_close(e1 * e1, one, "e1 e1 == 1");
    check_mv_close(e2 * e2, one, "e2 e2 == 1");
    check_mv_close(e12 * e12, -one, "e12 e12 == -1");
    check_mv_close(e1 * e2, e12, "e1 e2 == e12");
    check_mv_close(e2 * e1, -e12, "e2 e1 == -e12");

    // The pseudoscalar anticommutes with vectors. In Cl(3,0), e123 commutes
    // with everything; this is the line that makes Cl(2,0) its own algebra.
    check_mv_close(e1 * e12, e2, "e1 e12 == +e2");
    check_mv_close(e12 * e1, -e2, "e12 e1 == -e2");
    check_mv_close(e2 * e12, -e1, "e2 e12 == -e1");
    check_mv_close(e12 * e2, e1, "e12 e2 == +e1");
}

void test_cayley_table_against_permutation_rule()
{
    section("cayley table vs permutation rule");

    // All 16 basis products, each checked against the reference rule rather
    // than against a remembered table. This pins every sign in
    // geometric_product(), so a hand-edited term fails the build.
    int nonzero = 0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            int k = -1;
            const int sign = reference_blade_product(i, j, k);
            const Multivector<double> product = basis_blade(i) * basis_blade(j);
            const std::string what = std::string(kBladeNames[i]) + " * " + kBladeNames[j];
            check(k >= 0, what + " lands on a basis blade");
            if (k < 0) continue;
            ++nonzero;
            check_mv_close(product, static_cast<double>(sign) * basis_blade(k),
                           what + " == " + (sign > 0 ? "+" : "-") + kBladeNames[k]);
        }
    }
    // The metric is non-degenerate, so unlike Cl(3,0,1) no product vanishes.
    check(nonzero == 16, "all 16 basis products are non-zero");
}

void test_multivector_product()
{
    section("multivector product");

    const Multivector<double> a = make_mv(1.0, -2.0, 3.0, 0.5);
    const Multivector<double> b = make_mv(-0.5, 1.0, 2.0, -3.0);
    const Multivector<double> c = make_mv(2.0, 0.5, -1.0, 1.5);
    const Multivector<double> one = basis_blade(0);

    check_mv_close((a * b) * c, a * (b * c), "associativity", 1e-12);
    check_mv_close(a * (b + c), a * b + a * c, "left distributivity", 1e-12);
    check_mv_close((a + b) * c, a * c + b * c, "right distributivity", 1e-12);
    check_mv_close(one * a, a, "1 * a == a");
    check_mv_close(a * one, a, "a * 1 == a");
    check(mv_difference(a * b, b * a) > 1e-6, "the product is not commutative");
    check_mv_close(a * b, ga::geometric_product(a, b), "operator* == geometric_product");

    // Vector * Vector has a direct implementation; it must agree with the table.
    const Vector<double> u(1.0, 2.0), v(3.0, -4.0);
    check_mv_close(u * v, as_mv(u) * as_mv(v), "Vector * Vector == general product");
    check_mv_close(ga::geometric_product(u, v), u * v, "geometric_product(Vector, Vector)");
    const Rotor<double> uv = ga::rotor_product(u, v);
    check_rotor(uv, (u * v).scalar.value, (u * v).bivector.xy, "rotor_product(u, v)");

    // Rotor * Rotor is complex multiplication, also directly implemented.
    const Rotor<double> r(Scalar<double>(0.6), Bivector<double>(-0.8));
    const Rotor<double> s(Scalar<double>(2.0), Bivector<double>(0.5));
    const Multivector<double> rs = ga::to_multivector(r) * ga::to_multivector(s);
    check_rotor(r * s, rs.scalar.value, rs.bivector.xy, "Rotor * Rotor == general product");
    check_close(rs.vector.x, 0.0, "Rotor * Rotor has no vector part (x)");
    check_close(rs.vector.y, 0.0, "Rotor * Rotor has no vector part (y)");
    check_rotor_close(ga::rotor_product(r, s), r * s, "rotor_product(Rotor, Rotor)");

    // The even subalgebra of Cl(2,0) is commutative -- it is the complex
    // numbers. In Cl(3,0) rotors do NOT commute.
    check_rotor_close(r * s, s * r, "rotors commute in two dimensions");
}

void test_operator_return_types()
{
    section("operator* return types");

    using S = Scalar<double>;
    using V = Vector<double>;
    using B = Bivector<double>;
    using M = Multivector<double>;
    using R = Rotor<double>;

    // Every one of the 25 pairs, cell by cell as in docs/cl2.md. The narrow
    // cells are the point of the algebra and must not silently widen.
    static_assert(std::is_same<product_t<S, S>, S>::value, "S*S");
    static_assert(std::is_same<product_t<S, V>, V>::value, "S*V");
    static_assert(std::is_same<product_t<S, B>, B>::value, "S*B");
    static_assert(std::is_same<product_t<S, M>, M>::value, "S*M");
    static_assert(std::is_same<product_t<S, R>, R>::value, "S*R");

    static_assert(std::is_same<product_t<V, S>, V>::value, "V*S");
    static_assert(std::is_same<product_t<V, V>, M>::value, "V*V is a Multivector, as in Cl(3,0)");
    static_assert(std::is_same<product_t<V, B>, V>::value, "V*B is a pure Vector");
    static_assert(std::is_same<product_t<V, M>, M>::value, "V*M");
    static_assert(std::is_same<product_t<V, R>, V>::value, "V*R is a pure Vector");

    static_assert(std::is_same<product_t<B, S>, B>::value, "B*S");
    static_assert(std::is_same<product_t<B, V>, V>::value, "B*V is a pure Vector");
    static_assert(std::is_same<product_t<B, B>, S>::value, "B*B is a pure Scalar");
    static_assert(std::is_same<product_t<B, M>, M>::value, "B*M");
    static_assert(std::is_same<product_t<B, R>, R>::value, "B*R is a Rotor");

    static_assert(std::is_same<product_t<M, S>, M>::value, "M*S");
    static_assert(std::is_same<product_t<M, V>, M>::value, "M*V");
    static_assert(std::is_same<product_t<M, B>, M>::value, "M*B");
    static_assert(std::is_same<product_t<M, M>, M>::value, "M*M");
    static_assert(std::is_same<product_t<M, R>, M>::value, "M*R");

    static_assert(std::is_same<product_t<R, S>, R>::value, "R*S");
    static_assert(std::is_same<product_t<R, V>, V>::value, "R*V is a pure Vector");
    static_assert(std::is_same<product_t<R, B>, R>::value, "R*B is a Rotor");
    static_assert(std::is_same<product_t<R, M>, M>::value, "R*M");
    static_assert(std::is_same<product_t<R, R>, R>::value, "R*R");

    // Vector ^ Bivector would be grade 3, which does not exist here. It is
    // deliberately undefined rather than returning a type that pretends.
    static_assert(has_wedge<V, V>::value, "Vector ^ Vector is defined");
    static_assert(!has_wedge<V, B>::value, "Vector ^ Bivector must NOT be defined");
    static_assert(!has_wedge<B, V>::value, "Bivector ^ Vector must NOT be defined");

    check(true, "all 25 operator* return types and the wedge gap hold");
}

void test_mixed_products()
{
    section("mixed products");

    const Vector<double> v(3.0, -2.0);
    const Bivector<double> b(1.5);
    const Rotor<double> r(Scalar<double>(0.8), Bivector<double>(-0.6));
    const Multivector<double> m = make_mv(0.5, -1.0, 2.0, 0.25);

    const Multivector<double> vb = as_mv(v) * make_mv(0, 0, 0, b.xy);
    const Multivector<double> bv = make_mv(0, 0, 0, b.xy) * as_mv(v);
    const Multivector<double> vr = as_mv(v) * ga::to_multivector(r);
    const Multivector<double> rv = ga::to_multivector(r) * as_mv(v);
    const Multivector<double> br = make_mv(0, 0, 0, b.xy) * ga::to_multivector(r);
    const Multivector<double> rb = ga::to_multivector(r) * make_mv(0, 0, 0, b.xy);

    // Each narrow product equals the general product, AND the grades it drops
    // are zero there -- the narrowing loses nothing.
    check_vector_close(v * b, vb.vector, "Vector * Bivector");
    check_close(vb.scalar.value + vb.bivector.xy, 0.0, "Vector * Bivector drops nothing");
    check_vector_close(b * v, bv.vector, "Bivector * Vector");
    check_close(bv.scalar.value + bv.bivector.xy, 0.0, "Bivector * Vector drops nothing");
    check_vector_close(v * r, vr.vector, "Vector * Rotor");
    check_close(vr.scalar.value + vr.bivector.xy, 0.0, "Vector * Rotor drops nothing");
    check_vector_close(r * v, rv.vector, "Rotor * Vector");
    check_close(rv.scalar.value + rv.bivector.xy, 0.0, "Rotor * Vector drops nothing");
    check_rotor(b * r, br.scalar.value, br.bivector.xy, "Bivector * Rotor");
    check_close(br.vector.x + br.vector.y, 0.0, "Bivector * Rotor drops nothing");
    check_rotor(r * b, rb.scalar.value, rb.bivector.xy, "Rotor * Bivector");
    check_close(rb.vector.x + rb.vector.y, 0.0, "Rotor * Bivector drops nothing");

    // The pseudoscalar anticommutes with vectors but commutes with rotors.
    check_vector_close(b * v, -(v * b), "Bivector * Vector == -(Vector * Bivector)");
    check_rotor_close(b * r, r * b, "Bivector * Rotor == Rotor * Bivector");

    // Everything with a Multivector goes straight through the table.
    check_mv_close(m * v, m * as_mv(v), "Multivector * Vector");
    check_mv_close(v * m, as_mv(v) * m, "Vector * Multivector");
    check_mv_close(m * b, m * make_mv(0, 0, 0, b.xy), "Multivector * Bivector");
    check_mv_close(b * m, make_mv(0, 0, 0, b.xy) * m, "Bivector * Multivector");
    check_mv_close(m * r, m * ga::to_multivector(r), "Multivector * Rotor");
    check_mv_close(r * m, ga::to_multivector(r) * m, "Rotor * Multivector");

    // Bivector * Bivector (a member) agrees with the table's scalar part.
    const Bivector<double> b2(-4.0);
    check_scalar(b * b2, (make_mv(0, 0, 0, b.xy) * make_mv(0, 0, 0, b2.xy)).scalar.value,
                 "Bivector * Bivector == general product");
}

void test_dot_and_wedge()
{
    section("dot and wedge");

    const Vector<double> a(1.0, 2.0), b(3.0, 4.0);
    check_scalar(a | b, 11.0, "(1,2) | (3,4)");
    check_scalar(ga::dot_product(a, b), 11.0, "dot_product");
    check_bivector(a ^ b, -2.0, "(1,2) ^ (3,4)");
    check_bivector(ga::wedge_product(a, b), -2.0, "wedge_product");
    check_bivector(b ^ a, 2.0, "wedge is antisymmetric");
    check_bivector(a ^ a, 0.0, "v ^ v == 0");
    check_bivector(Vector<double>(1, 0) ^ Vector<double>(0, 1), 1.0, "e1 ^ e2 == +e12");

    // The geometric product of two vectors is exactly dot + wedge.
    const Multivector<double> ab = a * b;
    check_mv(ab, 11.0, 0.0, 0.0, -2.0, "a * b == (a|b) + (a^b)");
}

// ---------------------------------------------------------------------------
// Sums
// ---------------------------------------------------------------------------

void test_mixed_grade_addition()
{
    section("mixed-grade addition");

    const Scalar<double> s(1.0);
    const Vector<double> v(2.0, 3.0);
    const Bivector<double> b(4.0);
    const Multivector<double> m = make_mv(10.0, 20.0, 30.0, 40.0);

    check_mv(s + v, 1.0, 2.0, 3.0, 0.0, "Scalar + Vector");
    check_mv(v + s, 1.0, 2.0, 3.0, 0.0, "Vector + Scalar");
    check_mv(s + b, 1.0, 0.0, 0.0, 4.0, "Scalar + Bivector");
    check_mv(b + s, 1.0, 0.0, 0.0, 4.0, "Bivector + Scalar");
    check_mv(v + b, 0.0, 2.0, 3.0, 4.0, "Vector + Bivector");
    check_mv(b + v, 0.0, 2.0, 3.0, 4.0, "Bivector + Vector");
    check_mv(m + s, 11.0, 20.0, 30.0, 40.0, "Multivector + Scalar");
    check_mv(s + m, 11.0, 20.0, 30.0, 40.0, "Scalar + Multivector");
    check_mv(m + v, 10.0, 22.0, 33.0, 40.0, "Multivector + Vector");
    check_mv(v + m, 10.0, 22.0, 33.0, 40.0, "Vector + Multivector");
    check_mv(m + b, 10.0, 20.0, 30.0, 44.0, "Multivector + Bivector");
    check_mv(b + m, 10.0, 20.0, 30.0, 44.0, "Bivector + Multivector");
    check_mv(s + v + b, 1.0, 2.0, 3.0, 4.0, "chained sum");

    using S = Scalar<double>;
    using V = Vector<double>;
    using B = Bivector<double>;
    using M = Multivector<double>;
    using R = Rotor<double>;

    // Scalar + Bivector is exactly a rotor, but it widens anyway: one uniform
    // rule for every mixed sum.
    static_assert(std::is_same<sum_t<S, B>, M>::value, "Scalar + Bivector widens to Multivector");
    static_assert(std::is_same<sum_t<V, B>, M>::value, "Vector + Bivector");
    static_assert(std::is_same<sum_t<M, V>, M>::value, "Multivector + Vector");
    static_assert(std::is_same<sum_t<B, B>, B>::value, "same-grade sums keep their type");
    static_assert(std::is_same<sum_t<M, R>, M>::value, "Multivector + Rotor widens losslessly");

    // The named rotor spellings.
    const Rotor<double> r = ga::rotor_sum(s, b);
    check_rotor(r, 1.0, 4.0, "rotor_sum");
    check_rotor(ga::to_rotor(m), 10.0, 40.0, "to_rotor drops the vector part");
    check_mv(ga::to_multivector(r), 1.0, 0.0, 0.0, 4.0, "to_multivector");
    check_mv(m + r, 11.0, 20.0, 30.0, 44.0, "Multivector + Rotor");
}

void test_mixed_grade_subtraction()
{
    section("mixed-grade subtraction");

    const Scalar<double> s(1.0);
    const Vector<double> v(2.0, 3.0);
    const Bivector<double> b(4.0);
    const Multivector<double> m = make_mv(10.0, 20.0, 30.0, 40.0);

    check_mv(s - v, 1.0, -2.0, -3.0, 0.0, "Scalar - Vector");
    check_mv(v - s, -1.0, 2.0, 3.0, 0.0, "Vector - Scalar");
    check_mv(s - b, 1.0, 0.0, 0.0, -4.0, "Scalar - Bivector");
    check_mv(b - s, -1.0, 0.0, 0.0, 4.0, "Bivector - Scalar");
    check_mv(v - b, 0.0, 2.0, 3.0, -4.0, "Vector - Bivector");
    check_mv(b - v, 0.0, -2.0, -3.0, 4.0, "Bivector - Vector");
    check_mv(m - s, 9.0, 20.0, 30.0, 40.0, "Multivector - Scalar");
    check_mv(s - m, -9.0, -20.0, -30.0, -40.0, "Scalar - Multivector");
    check_mv(m - v, 10.0, 18.0, 27.0, 40.0, "Multivector - Vector");
    check_mv(v - m, -10.0, -18.0, -27.0, -40.0, "Vector - Multivector");
    check_mv(m - b, 10.0, 20.0, 30.0, 36.0, "Multivector - Bivector");
    check_mv(b - m, -10.0, -20.0, -30.0, -36.0, "Bivector - Multivector");

    static_assert(std::is_same<difference_t<Scalar<double>, Vector<double>>, Multivector<double>>::value,
                  "Scalar - Vector widens to Multivector");
    static_assert(std::is_same<difference_t<Vector<double>, Vector<double>>, Vector<double>>::value,
                  "same-grade differences keep their type");
}

void test_grade_projection()
{
    section("grade projection");

    const Multivector<double> m = make_mv(1.0, 2.0, 3.0, 4.0);
    check_scalar(ga::grade0(m), 1.0, "grade0(Multivector)");
    check_vector(ga::grade1(m), 2.0, 3.0, "grade1(Multivector)");
    check_bivector(ga::grade2(m), 4.0, "grade2(Multivector)");

    const Rotor<double> r(Scalar<double>(5.0), Bivector<double>(6.0));
    check_scalar(ga::grade0(r), 5.0, "grade0(Rotor)");
    check_bivector(ga::grade2(r), 6.0, "grade2(Rotor)");

    // The overload Cl(3,0) lacks: grade0 of a Scalar is the identity, which is
    // what lets scalar_product() accept two bivectors.
    check_scalar(ga::grade0(Scalar<double>(7.0)), 7.0, "grade0(Scalar)");
}

void test_contractions()
{
    section("contractions");

    const Vector<double> u(1.0, 2.0), v(3.0, -1.0);
    const Bivector<double> b(2.0), c(-3.0);

    // Each contraction is the matching grade of the general product.
    const Multivector<double> vb = as_mv(v) * make_mv(0, 0, 0, b.xy);
    const Multivector<double> bv = make_mv(0, 0, 0, b.xy) * as_mv(v);
    const Multivector<double> bc = make_mv(0, 0, 0, b.xy) * make_mv(0, 0, 0, c.xy);

    check_scalar(ga::left_contraction(u, v), (as_mv(u) * as_mv(v)).scalar.value, "left(V, V)");
    check_vector_close(ga::left_contraction(v, b), vb.vector, "left(V, B)");
    check_scalar(ga::left_contraction(b, c), bc.scalar.value, "left(B, B)");
    check_vector_close(v | b, vb.vector, "V | B");
    check_scalar(b | c, bc.scalar.value, "B | B");
    check_scalar(b | c, 6.0, "(2 e12) | (-3 e12) == 6");

    check_scalar(ga::right_contraction(u, v), (u | v).value, "right(V, V) == left(V, V)");
    check_vector_close(ga::right_contraction(b, v), bv.vector, "right(B, V)");
    check_scalar(ga::right_contraction(b, c), bc.scalar.value, "right(B, B)");

    // e12 anticommutes with vectors, so right and left differ by a sign.
    check_vector_close(ga::right_contraction(b, v), -ga::left_contraction(v, b),
                       "right(B, V) == -left(V, B)");

    // e1 contracted into e12 is e2: the quarter turn again.
    check_vector(Vector<double>(1, 0) | Bivector<double>(1), 0.0, 1.0, "e1 | e12 == e2");

    // scalar_product works for every pair whose product has a grade 0 part.
    const Rotor<double> r(Scalar<double>(0.6), Bivector<double>(0.8));
    const Rotor<double> s(Scalar<double>(2.0), Bivector<double>(-1.0));
    const Multivector<double> m = make_mv(1.0, 2.0, 3.0, 4.0);
    const Multivector<double> n = make_mv(-1.0, 0.5, 2.0, -3.0);
    check_scalar(ga::scalar_product(u, v), (u | v).value, "scalar_product(V, V)");
    check_scalar(ga::scalar_product(b, c), (b * c).value, "scalar_product(B, B)");
    check_scalar(ga::scalar_product(r, s), (r * s).scalar.value, "scalar_product(R, R)");
    check_scalar(ga::scalar_product(m, n), (m * n).scalar.value, "scalar_product(M, M)");
}

// ---------------------------------------------------------------------------
// Magnitudes and involutions
// ---------------------------------------------------------------------------

void test_norms()
{
    section("norms");

    check_scalar(ga::norm(Scalar<double>(-2.0)), 2.0, "norm(Scalar)");
    check_scalar(ga::squared_norm(Scalar<double>(-2.0)), 4.0, "squared_norm(Scalar)");
    check_scalar(ga::norm(Vector<double>(3.0, 4.0)), 5.0, "norm(Vector)");
    check_scalar(ga::squared_norm(Vector<double>(3.0, 4.0)), 25.0, "squared_norm(Vector)");
    check_scalar(ga::norm(Bivector<double>(-6.0)), 6.0, "norm(Bivector)");
    check_scalar(ga::squared_norm(Bivector<double>(-6.0)), 36.0, "squared_norm(Bivector)");
    check_scalar(ga::norm(Rotor<double>(Scalar<double>(3.0), Bivector<double>(4.0))), 5.0, "norm(Rotor)");
    check_scalar(ga::squared_norm(Rotor<double>(Scalar<double>(3.0), Bivector<double>(4.0))), 25.0,
                 "squared_norm(Rotor)");

    // squared_norm is the grade 0 part of m * reverse(m), and the metric is
    // Euclidean, so every component counts.
    const Multivector<double> m = make_mv(1.0, -2.0, 3.0, 0.5);
    check_scalar(ga::squared_norm(m), 1.0 + 4.0 + 9.0 + 0.25, "squared_norm(Multivector)");
    check_scalar(ga::squared_norm(m), (m * ga::reverse(m)).scalar.value, "squared_norm == <m ~m>_0");
    check_scalar(ga::norm(m), std::sqrt(14.25), "norm(Multivector)");
}

void test_normalize()
{
    section("normalize");

    check_scalar(ga::normalize(Scalar<double>(-3.0)), -1.0, "normalize(Scalar) keeps the sign");
    check_vector(ga::normalize(Vector<double>(3.0, 4.0)), 0.6, 0.8, "normalize(Vector)");
    // A bivector has one component, so its unit version is +e12 or -e12.
    check_bivector(ga::normalize(Bivector<double>(-5.0)), -1.0, "normalize(Bivector) keeps orientation");
    check_rotor(ga::normalize(Rotor<double>(Scalar<double>(3.0), Bivector<double>(4.0))), 0.6, 0.8,
                "normalize(Rotor)");
    check_close(ga::norm(ga::normalize(make_mv(1.0, -2.0, 3.0, 0.5))).value, 1.0,
                "normalize(Multivector) is unit");
}

void test_involutions()
{
    section("involutions");

    const Multivector<double> m = make_mv(1.0, 2.0, 3.0, 4.0);

    // The sign table in involutions.hpp, grade by grade.
    check_mv(ga::reverse(m), 1.0, 2.0, 3.0, -4.0, "reverse: + + -");
    check_mv(ga::involute(m), 1.0, -2.0, -3.0, 4.0, "involute: + - +");
    check_mv(ga::conjugate(m), 1.0, -2.0, -3.0, -4.0, "conjugate: + - -");

    // conjugate is reverse composed with involute, in either order.
    check_mv_close(ga::conjugate(m), ga::reverse(ga::involute(m)), "conjugate == reverse . involute");
    check_mv_close(ga::conjugate(m), ga::involute(ga::reverse(m)), "conjugate == involute . reverse");

    // The per-type overloads agree with the multivector one.
    check_scalar(ga::reverse(Scalar<double>(2.0)), 2.0, "reverse(Scalar)");
    check_vector(ga::reverse(Vector<double>(1.0, 2.0)), 1.0, 2.0, "reverse(Vector)");
    check_bivector(ga::reverse(Bivector<double>(3.0)), -3.0, "reverse(Bivector)");
    check_rotor(ga::reverse(Rotor<double>(Scalar<double>(1.0), Bivector<double>(2.0))), 1.0, -2.0,
                "reverse(Rotor)");
    check_scalar(ga::involute(Scalar<double>(2.0)), 2.0, "involute(Scalar)");
    check_vector(ga::involute(Vector<double>(1.0, 2.0)), -1.0, -2.0, "involute(Vector)");
    check_bivector(ga::involute(Bivector<double>(3.0)), 3.0, "involute(Bivector)");
    check_rotor(ga::involute(Rotor<double>(Scalar<double>(1.0), Bivector<double>(2.0))), 1.0, 2.0,
                "involute(Rotor)");
    check_scalar(ga::conjugate(Scalar<double>(2.0)), 2.0, "conjugate(Scalar)");
    check_vector(ga::conjugate(Vector<double>(1.0, 2.0)), -1.0, -2.0, "conjugate(Vector)");
    check_bivector(ga::conjugate(Bivector<double>(3.0)), -3.0, "conjugate(Bivector)");
    check_rotor(ga::conjugate(Rotor<double>(Scalar<double>(1.0), Bivector<double>(2.0))), 1.0, -2.0,
                "conjugate(Rotor) == reverse(Rotor)");

    // The anti-automorphism and automorphism properties, on full multivectors.
    const Multivector<double> n = make_mv(-0.5, 1.5, -2.0, 0.75);
    check_mv_close(ga::reverse(m * n), ga::reverse(n) * ga::reverse(m), "reverse(ab) == ~b ~a");
    check_mv_close(ga::involute(m * n), ga::involute(m) * ga::involute(n),
                   "involute(ab) == involute(a) involute(b)");
    check_mv_close(ga::conjugate(m * n), ga::conjugate(n) * ga::conjugate(m),
                   "conjugate(ab) == conj(b) conj(a)");
}

void test_inverse()
{
    section("inverse");

    const Multivector<double> one = basis_blade(0);

    check_scalar(Scalar<double>(4.0) * ga::inverse(Scalar<double>(4.0)), 1.0, "s * inverse(s)");

    const Vector<double> v(3.0, -4.0);
    check_scalar((v * ga::inverse(v)).scalar, 1.0, "v * inverse(v)");
    check_vector(ga::inverse(Vector<double>(1, 0)), 1.0, 0.0, "inverse(e1) == e1");

    // e12^2 = -1, so the bivector inverse carries a minus sign.
    const Bivector<double> b(2.5);
    check_scalar(b * ga::inverse(b), 1.0, "b * inverse(b)");
    check_bivector(ga::inverse(Bivector<double>(1.0)), -1.0, "inverse(e12) == -e12");

    // A non-unit rotor, so the division by |r|^2 is exercised.
    const Rotor<double> r(Scalar<double>(3.0), Bivector<double>(4.0));
    check_rotor(r * ga::inverse(r), 1.0, 0.0, "r * inverse(r)");
    check_rotor(ga::inverse(r) * r, 1.0, 0.0, "inverse(r) * r");

    // The fact the multivector inverse rests on: m * conjugate(m) is a pure
    // scalar in Cl(2,0), and it is a^2 - x^2 - y^2 + b^2.
    const Multivector<double> m = make_mv(1.0, -2.0, 3.0, 0.5);
    const Multivector<double> collapsed = m * ga::conjugate(m);
    check_mv(collapsed, 1.0 - 4.0 - 9.0 + 0.25, 0.0, 0.0, 0.0, "m * conjugate(m) is a pure scalar");

    check_mv_close(m * ga::inverse(m), one, "m * inverse(m) == 1");
    check_mv_close(ga::inverse(m) * m, one, "inverse(m) * m == 1");

    // That denominator is indefinite. Here it is positive...
    const Multivector<double> p = make_mv(3.0, 1.0, 1.0, 2.0);
    check_mv_close(p * ga::inverse(p), one, "invertible with a positive denominator");
    check_mv_close(ga::inverse(p) * p, one, "invertible with a positive denominator (left)");

    // ...and a non-zero multivector can still be singular. (1 + e1)(1 - e1) = 0.
    const Multivector<double> singular = make_mv(1.0, 1.0, 0.0, 0.0);
    check_mv_close(singular * ga::conjugate(singular), Multivector<double>(),
                   "(1 + e1) is a zero divisor");
    check_mv_close(ga::inverse(singular), Multivector<double>(), "inverse(1 + e1) is zero");
}

void test_dual()
{
    section("dual");

    // The table in dual.hpp: 1 -> e12, e1 -> e2, e2 -> -e1, e12 -> -1.
    check_bivector(ga::dual(Scalar<double>(1.0)), 1.0, "dual(1) == e12");
    check_vector(ga::dual(Vector<double>(1, 0)), 0.0, 1.0, "dual(e1) == e2");
    check_vector(ga::dual(Vector<double>(0, 1)), -1.0, 0.0, "dual(e2) == -e1");
    check_scalar(ga::dual(Bivector<double>(1.0)), -1.0, "dual(e12) == -1");

    // dual(A) is A * e12, checked against the general product.
    const Bivector<double> I(1.0);
    const Vector<double> v(3.0, -2.0);
    const Multivector<double> m = make_mv(1.0, 2.0, 3.0, 4.0);
    check_vector_close(ga::dual(v), v * I, "dual(v) == v * e12");
    check_mv_close(ga::dual(m), m * I, "dual(m) == m * e12");

    // Applying it twice negates, on every grade. Correct, not a defect.
    check_scalar(ga::dual(ga::dual(Scalar<double>(5.0))), -5.0, "dual(dual(s)) == -s");
    check_vector(ga::dual(ga::dual(v)), -3.0, 2.0, "dual(dual(v)) == -v");
    check_bivector(ga::dual(ga::dual(Bivector<double>(5.0))), -5.0, "dual(dual(b)) == -b");
    check_mv_close(ga::dual(ga::dual(m)), -m, "dual(dual(m)) == -m");
}

void test_comparison()
{
    section("comparison");

    // The default tolerance is 100 epsilons, about 2.2e-14 for double: enough
    // to absorb rounding, and deliberately too small to hide a 1e-10 error.
    const Vector<double> a(1.0, 2.0), b(1.0, 2.0);
    const Vector<double> rounding(1.0, 2.0 + 1e-15), error(1.0, 2.0 + 1e-10);
    check(a == b, "exact == on equal vectors");
    check(a != rounding, "exact != on a 1e-15 difference");
    check(ga::approx_equal(a, rounding), "approx_equal absorbs rounding by default");
    check(!ga::approx_equal(a, error), "approx_equal does not absorb a 1e-10 error by default");
    check(ga::approx_equal(a, error, 1e-9), "approx_equal with an explicit, looser tolerance");

    check(Scalar<double>(1.0) == Scalar<double>(1.0), "Scalar ==");
    check(Scalar<double>(1.0) != Scalar<double>(2.0), "Scalar !=");
    check(Bivector<double>(3.0) == Bivector<double>(3.0), "Bivector ==");
    check(Bivector<double>(3.0) != Bivector<double>(-3.0), "Bivector !=");
    check(make_mv(1, 2, 3, 4) == make_mv(1, 2, 3, 4), "Multivector ==");
    check(make_mv(1, 2, 3, 4) != make_mv(1, 2, 3, 5), "Multivector !=");
    check(ga::approx_equal(make_mv(1, 2, 3, 4), make_mv(1, 2, 3, 4 + 1e-15)), "Multivector approx_equal");
    check(ga::approx_equal(Scalar<double>(1.0), Scalar<double>(1.0 + 1e-15)), "Scalar approx_equal");
    check(ga::approx_equal(Bivector<double>(1.0), Bivector<double>(1.0 + 1e-15)), "Bivector approx_equal");

    // A rotor and its negation are the same rotation, but not the same rotor.
    const Rotor<double> r(Scalar<double>(0.6), Bivector<double>(0.8));
    const Rotor<double> neg(Scalar<double>(-0.6), Bivector<double>(-0.8));
    check(r == r, "Rotor ==");
    check(r != neg, "a rotor is not equal to its negation");
    check(!ga::approx_equal(r, neg), "nor approx_equal to it");
    check(ga::approx_equal(r, Rotor<double>(Scalar<double>(0.6 + 1e-15), Bivector<double>(0.8))),
          "Rotor approx_equal");

    // operator<< delegates to to_string() for every type.
    std::ostringstream os;
    os << Scalar<double>(1.5) << '|' << Vector<double>(1, 2) << '|' << Bivector<double>(3)
       << '|' << r << '|' << make_mv(1, 2, 3, 4);
    check(os.str() == "1.5|1*e1 + 2*e2|3*e12|" + r.to_string() + "|(1) + (2*e1 + 3*e2) + (4*e12)",
          "operator<< matches to_string for every type");
}

// ---------------------------------------------------------------------------
// Rotations
// ---------------------------------------------------------------------------

void test_exp_log()
{
    section("exp and log");

    check_rotor(ga::exp(Bivector<double>(0.0)), 1.0, 0.0, "exp(0) is the identity");

    // Euler's formula, with e12 as i.
    for (int i = -30; i <= 30; ++i) {
        const double theta = 0.1 * i;
        const Rotor<double> r = ga::exp(Bivector<double>(theta));
        const std::string at = " at theta = " + std::to_string(theta);
        check_rotor(r, std::cos(theta), std::sin(theta), "exp(theta e12) == cos + sin e12" + at);
        check_close(ga::norm(r).value, 1.0, "exp is always unit" + at);
        check_bivector(ga::log(r), theta, "log(exp(b)) == b" + at);
    }

    // log undoes exp, and exp undoes log, on a unit rotor.
    const Rotor<double> r = ga::normalize(Rotor<double>(Scalar<double>(-0.3), Bivector<double>(0.9)));
    check_rotor_close(ga::exp(ga::log(r)), r, "exp(log(r)) == r");

    // A negative scalar part means a turn past a half: acos-based logs need a
    // clamp and a sign fix here, atan2 does not.
    check_bivector(ga::log(ga::exp(Bivector<double>(2.5))), 2.5, "log of a turn past a half");

    // atan2 only sees the ratio, so drift off unit length does not move the angle.
    const Rotor<double> drifted = r * Scalar<double>(3.7);
    check_bivector(ga::log(drifted), ga::log(r).xy, "log ignores a rotor's magnitude");
    check_bivector(ga::log(Rotor<double>()), 0.0, "log of the zero rotor is zero");

    // The range is (-pi, pi]: a -0.0 bivector part must not push -1 to -pi.
    check_bivector(ga::log(Rotor<double>(Scalar<double>(-1.0), Bivector<double>(-0.0))), kPi,
                   "log(-1 - 0.0 e12) == +pi e12");

    // exp(theta B) rotates by 2 theta -- clockwise for positive theta, which is
    // why rotor_from_angle carries a minus sign.
    const Vector<double> turned = ga::sandwich(Vector<double>(1, 0), ga::exp(Bivector<double>(0.4)));
    check_close(std::atan2(turned.y, turned.x), -0.8, "exp(0.4 e12) turns e1 by -0.8");
}

void test_rotor_construction()
{
    section("rotor construction");

    const Vector<double> e1(1, 0), e2(0, 1);

    check_rotor(ga::identity_rotor<double>(), 1.0, 0.0, "identity_rotor");
    check_vector_close(ga::rotate(Vector<double>(3, -2), ga::identity_rotor<double>()),
                       Vector<double>(3, -2), "the identity rotates nothing");

    // Counter-clockwise for a positive angle, the same sign as Cl(3,0) about +z.
    const Rotor<double> quarter = ga::rotor_from_angle(kPi / 2);
    check_rotor(quarter, std::cos(kPi / 4), -std::sin(kPi / 4), "rotor_from_angle(pi/2) components");
    check_vector_close(ga::rotate(e1, quarter), e2, "+90 degrees takes e1 to e2");
    check_vector_close(ga::rotate(e2, quarter), -e1, "+90 degrees takes e2 to -e1");
    check_vector_close(ga::rotate(e1, ga::rotor_from_angle(-kPi / 2)), -e2, "-90 degrees takes e1 to -e2");
    check_vector_close(ga::rotate(e1, ga::rotor_from_angle(kPi)), -e1, "180 degrees takes e1 to -e1");

    // Rotations add, and in two dimensions they commute.
    check_rotor_close(ga::rotor_from_angle(0.3) * ga::rotor_from_angle(0.9), ga::rotor_from_angle(1.2),
                      "rotor_from_angle(a) * rotor_from_angle(b) == rotor_from_angle(a + b)");

    // rotor_angle reads the angle back, in (-pi, pi].
    for (int i = -30; i <= 30; ++i) {
        const double theta = 0.1 * i;
        check_close(ga::rotor_angle(ga::rotor_from_angle(theta)), theta,
                    "rotor_angle round trip at " + std::to_string(theta));
    }
    const Rotor<double> r = ga::rotor_from_angle(1.1);
    check_close(ga::rotor_angle(r * Scalar<double>(-1.0)), 1.1, "rotor_angle(-r) == rotor_angle(r)");
    check_close(ga::rotor_angle(r * Scalar<double>(2.5)), 1.1, "rotor_angle ignores magnitude");
    check_close(ga::rotor_angle(ga::rotor_from_angle(1.5 * kPi)), -0.5 * kPi,
                "rotor_angle reports the equivalent angle in range");
    check_close(ga::rotor_angle(ga::identity_rotor<double>()), 0.0, "rotor_angle(identity) == 0");

    // Both exact half-turn rotors report +pi. For 0 + 1 e12 the numerator is
    // (-2 * 0) * 1 = -0.0, which plain atan2 would turn into -pi.
    check_close(ga::rotor_angle(Rotor<double>(Scalar<double>(0), Bivector<double>(1))), kPi,
                "rotor_angle(+e12) == +pi");
    check_close(ga::rotor_angle(Rotor<double>(Scalar<double>(0), Bivector<double>(-1))), kPi,
                "rotor_angle(-e12) == +pi");

    // rotor_between lands exactly, from non-unit inputs.
    const Vector<double> from(3.0, 1.0), to(-1.0, 2.0);
    const Rotor<double> between = ga::rotor_between(from, to);
    check_close(ga::norm(between).value, 1.0, "rotor_between is unit");
    check_vector_close(ga::rotate(ga::normalize(from), between), ga::normalize(to),
                       "rotor_between carries from onto to");
    check_close(ga::rotor_angle(ga::rotor_between(Vector<double>(1, 0), Vector<double>(0, 2))), kPi / 2,
                "rotor_between(e1, 2 e2) is a quarter turn");

    // Opposite: no plane to choose in 2D, only a direction. It is the exact
    // counter-clockwise half turn.
    const Rotor<double> half = ga::rotor_between(from, -from);
    check_rotor(half, 0.0, -1.0, "rotor_between of opposites is rotor_from_angle(pi)");
    check_vector_close(ga::rotate(from, half), -from, "and it lands on the opposite");

    check_rotor(ga::rotor_between(Vector<double>(), e1), 1.0, 0.0, "a zero input gives the identity");
}

void test_sandwich()
{
    section("sandwich");

    const Vector<double> v(2.0, -1.5);
    const double scales[3] = {0.5, 1.0, 2.5};

    for (int i = -30; i <= 30; ++i) {
        for (double k : scales) {
            // Non-unit rotors too: every identity below holds at any magnitude.
            const Rotor<double> r = ga::rotor_from_angle(0.1 * i) * Scalar<double>(k);
            const std::string at = " (angle " + std::to_string(0.1 * i) + ", scale "
                                 + std::to_string(k) + ")";

            // The closed form agrees with the two-sided general product.
            const Multivector<double> general =
                ga::to_multivector(r) * as_mv(v) * ga::to_multivector(ga::reverse(r));
            check_vector_close(ga::sandwich(v, r), general.vector, "closed form == R v ~R" + at, 1e-11);
            check_close(general.scalar.value, 0.0, "R v ~R has no scalar part" + at, 1e-11);
            check_close(general.bivector.xy, 0.0, "R v ~R has no bivector part" + at, 1e-11);

            // The 2D-only shortcut: the sandwich is a one-sided product.
            check_vector_close(ga::sandwich(v, r), r * r * v, "R v ~R == R^2 v" + at, 1e-11);
        }
    }

    // Rotation preserves length and the angle between two vectors.
    const Rotor<double> r = ga::rotor_from_angle(0.7);
    const Vector<double> w(-0.5, 3.0);
    check_close(ga::norm(ga::rotate(v, r)).value, ga::norm(v).value, "rotation preserves length");
    check_close(ga::angle_between(ga::rotate(v, r), ga::rotate(w, r)), ga::angle_between(v, w),
                "rotation preserves angles");
    check_close(ga::angle_between(v, ga::rotate(v, r)), 0.7, "rotate turns by the rotor's angle");

    // Composition: rotating by r1 then r2 is rotating by r2 * r1.
    const Rotor<double> r1 = ga::rotor_from_angle(0.4), r2 = ga::rotor_from_angle(-1.3);
    check_vector_close(ga::rotate(ga::rotate(v, r1), r2), ga::rotate(v, r2 * r1),
                       "rotate(rotate(v, r1), r2) == rotate(v, r2 r1)");

    // There is only one plane, so a unit rotor leaves a bivector alone, and a
    // rotor of norm k scales it by k^2.
    const Bivector<double> b(3.0);
    check_bivector(ga::sandwich(b, r), 3.0, "sandwich(B, unit r) == B");
    check_bivector(ga::rotate(b, r), 3.0, "rotate(B, r) == B");
    check_bivector(ga::sandwich(b, r * Scalar<double>(2.0)), 12.0, "sandwich(B, 2r) == 4B");

    // A multivector rotates grade by grade.
    const Multivector<double> m = make_mv(1.5, v.x, v.y, -2.0);
    const Multivector<double> rm = ga::sandwich(m, r);
    check_scalar(rm.scalar, 1.5, "sandwich(m, r) keeps the scalar");
    check_vector_close(rm.vector, ga::sandwich(v, r), "sandwich(m, r) rotates the vector");
    check_bivector(rm.bivector, -2.0, "sandwich(m, r) keeps the bivector");
    check_mv_close(ga::rotate(m, r), rm, "rotate(m, r) == sandwich(m, r)");
}

void test_slerp()
{
    section("slerp");

    const Rotor<double> a = ga::rotor_from_angle(0.2);
    const Rotor<double> b = ga::rotor_from_angle(1.0);

    check_rotor_close(ga::slerp(a, b, 0.0), a, "slerp at t = 0");
    check_rotor_close(ga::slerp(a, b, 1.0), b, "slerp at t = 1");
    check_close(ga::rotor_angle(ga::slerp(a, b, 0.25)), 0.4, "slerp interpolates the angle");
    check_close(ga::rotor_angle(ga::slerp(a, b, 0.5)), 0.6, "slerp midpoint");

    // -b is the same rotation as b; slerp must not take the long way round.
    check_close(ga::rotor_angle(ga::slerp(a, b * Scalar<double>(-1.0), 0.25)), 0.4,
                "slerp to -b takes the short arc");

    // Across the +-pi seam: from -3.0 to +3.0 the short way is through pi,
    // a turn of 2 pi - 6, not the 6 radians the other way.
    const Rotor<double> mid = ga::slerp(ga::rotor_from_angle(-3.0), ga::rotor_from_angle(3.0), 0.5);
    check_close(std::cos(ga::rotor_angle(mid)), -1.0, "slerp across the seam goes through pi", 1e-12);
    check_close(ga::norm(mid).value, 1.0, "slerp result is unit");
}

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

void test_geometry()
{
    section("geometry");

    const Vector<double> e1(1, 0), e2(0, 1);
    const Vector<double> v(3.0, 1.0);
    const Vector<double> diagonal(1.0, 1.0);

    // reflect(v, n) mirrors in the line perpendicular to n. Values from the
    // independent model.
    check_vector(ga::reflect(e1, e1), -1.0, 0.0, "reflect(e1, n = e1)");
    check_vector(ga::reflect(e2, e1), 0.0, 1.0, "reflect(e2, n = e1)");
    check_vector(ga::reflect(v, e1), -3.0, 1.0, "reflect((3,1), n = e1)");
    check_vector(ga::reflect(v, e2), 3.0, -1.0, "reflect((3,1), n = e2)");
    check_vector(ga::reflect(e1, diagonal), 0.0, -1.0, "reflect(e1, n = (1,1))");
    check_vector(ga::reflect(v, diagonal), -1.0, -3.0, "reflect((3,1), n = (1,1))");
    check_vector_close(ga::reflect(v, diagonal * Scalar<double>(5.0)), ga::reflect(v, diagonal),
                       "the normal need not be unit");
    check_vector_close(ga::reflect(ga::reflect(v, diagonal), diagonal), v, "reflecting twice is the identity");

    // Reflecting in the line ALONG u is -reflect(v, u).
    check_vector(-ga::reflect(v, e1), 3.0, -1.0, "-reflect(v, e1) mirrors in the x axis");

    // Two reflections make a rotation, by twice the angle from the first
    // normal to the second.
    const Vector<double> n2(std::cos(0.3), std::sin(0.3));
    check_vector_close(ga::reflect(ga::reflect(v, e1), n2), ga::rotate(v, ga::rotor_from_angle(0.6)),
                       "two reflections == a rotation by twice the angle");

    // Projection and rejection split a vector exactly.
    const Vector<double> u(2.0, 0.0);
    check_vector(ga::project(Vector<double>(3, 4), u), 3.0, 0.0, "project onto the x axis");
    check_vector(ga::reject(Vector<double>(3, 4), u), 0.0, 4.0, "reject from the x axis");
    const Vector<double> t(1.0, -2.0);
    check_vector_close(ga::project(v, t) + ga::reject(v, t), v, "project + reject == v");
    check_close((ga::reject(v, t) | t).value, 0.0, "reject is perpendicular to u");
    check_close((ga::project(v, t) ^ t).xy, 0.0, "project is parallel to u");

    // perp is the counter-clockwise quarter turn, and it is exactly dual().
    check_vector(ga::perp(Vector<double>(3, 4)), -4.0, 3.0, "perp((3,4)) == (-4,3)");
    check_vector_close(ga::perp(v), ga::dual(v), "perp == dual");
    check_close((ga::perp(v) | v).value, 0.0, "perp(v) is perpendicular to v");
    check_vector_close(ga::perp(ga::perp(v)), -v, "perp twice negates");
    check_vector_close(ga::perp(v), ga::rotate(v, ga::rotor_from_angle(kPi / 2)), "perp is +90 degrees");

    // signed_area is positive counter-clockwise, and is the wedge component.
    check_scalar(ga::signed_area(e1, e2), 1.0, "signed_area(e1, e2) == +1");
    check_scalar(ga::signed_area(e2, e1), -1.0, "signed_area(e2, e1) == -1");
    check_scalar(ga::signed_area(Vector<double>(1, 2), Vector<double>(3, 4)), -2.0, "signed_area((1,2),(3,4))");
    check_scalar(ga::signed_area(v, v), 0.0, "signed_area(v, v) == 0");
    check_scalar(ga::signed_area(v, ga::perp(v)), ga::squared_norm(v).value, "signed_area(v, perp v) == |v|^2");

    // angle_between is signed; angle is the direction from +x.
    check_close(ga::angle_between(e1, e2), kPi / 2, "angle_between(e1, e2) == +pi/2");
    check_close(ga::angle_between(e2, e1), -kPi / 2, "angle_between(e2, e1) == -pi/2");
    check_close(ga::angle_between(e1, -e1), kPi, "angle_between(e1, -e1) == pi");
    check_close(ga::angle_between(e1 * Scalar<double>(4.0), diagonal), kPi / 4, "angle_between ignores length");
    check_close(ga::angle(e1), 0.0, "angle(e1) == 0");
    check_close(ga::angle(e2), kPi / 2, "angle(e2) == pi/2");
    check_close(ga::angle(-e1), kPi, "angle(-e1) == pi");
    check_close(ga::angle(-e2), -kPi / 2, "angle(-e2) == -pi/2");

    // The principal range is (-pi, pi], and signed zero is what threatens it:
    // each input below makes plain atan2 return -pi. -e1 is (-1, -0.0), and
    // (-1, 0) ^ (1, 0) is -0.0.
    check(ga::angle(-e1) > 0.0, "angle(-e1) is +pi, not -pi");
    check(ga::angle_between(Vector<double>(-1, 0), e1) > 0.0, "a half turn between vectors is +pi");
    check_close(ga::angle_between(Vector<double>(-1, 0), e1), kPi, "angle_between((-1,0), e1) == pi");
}

// ---------------------------------------------------------------------------
// The even subalgebra is the complex numbers
// ---------------------------------------------------------------------------

void test_complex_isomorphism()
{
    section("complex isomorphism");

    // std::complex is an independent implementation of the even subalgebra, so
    // every rotor operation is checked against it. The alias makes the reading
    // explicit.
    using C = std::complex<double>;
    const auto as_complex = [](const Rotor<double>& r) { return C(r.scalar.value, r.bivector.xy); };

    const ga::Complex<double> z1(Scalar<double>(1.5), Bivector<double>(-0.5));
    const ga::Complex<double> z2(Scalar<double>(-0.25), Bivector<double>(2.0));
    const C c1 = as_complex(z1), c2 = as_complex(z2);

    const C product = c1 * c2;
    check_rotor(z1 * z2, product.real(), product.imag(), "rotor product == complex product");

    const C quotient = C(1.0, 0.0) / c1;
    check_rotor(ga::inverse(z1), quotient.real(), quotient.imag(), "inverse == 1/z");

    check_close(ga::norm(z1).value, std::abs(c1), "norm == |z|");
    check_rotor(ga::reverse(z1), std::conj(c1).real(), std::conj(c1).imag(), "reverse == conj(z)");
    check_close(ga::log(z1).xy, std::arg(c1), "log == arg(z)");

    const C e = std::exp(C(0.0, 0.7));
    check_rotor(ga::exp(Bivector<double>(0.7)), e.real(), e.imag(), "exp(theta e12) == e^(i theta)");

    // And a rotation in the plane is "multiply by e^(i theta)" -- here as R^2.
    const C point(2.0, -1.5);
    const C turned = std::exp(C(0.0, 1.1)) * point;
    const Vector<double> rotated = ga::rotate(Vector<double>(2.0, -1.5), ga::rotor_from_angle(1.1));
    check_close(rotated.x, turned.real(), "rotate == multiply by e^(i theta) (x)");
    check_close(rotated.y, turned.imag(), "rotate == multiply by e^(i theta) (y)");
}

// ---------------------------------------------------------------------------
// Instantiation sweep
// ---------------------------------------------------------------------------

// Templates only fail on instantiation, so including a header proves nothing.
// This calls every public entry point that exists so far, for each component
// type. It grows as operations land.
template<typename T>
void instantiate_every_entry_point()
{
    Scalar<T> s(T(2));
    Vector<T> v(T(1), T(2));
    Bivector<T> b(T(3));
    Rotor<T> r(s, b);
    Multivector<T> m(s, v, b);

    (void)(s + s); (void)(s - s); (void)(s * s); (void)(s / s); (void)(-s);
    (void)(v + v); (void)(v - v); (void)(-v); (void)v.magnitude();
    (void)(b + b); (void)(b - b); (void)(-b); (void)b.magnitude();
    (void)(b * b); (void)(b / b); (void)(b * T(2)); (void)(b / T(2)); (void)(T(2) * b);
    (void)(r + r); (void)(r - r);
    (void)(m + m); (void)(m - m); (void)(-m);

    (void)(s * v); (void)(v * s); (void)(v / s); (void)(T(2) * v);
    (void)(s * b); (void)(b * s); (void)(b / s);
    (void)(s * r); (void)(r * s); (void)(r / s); (void)(T(2) * r);
    (void)(s * m); (void)(m * s); (void)(m / s); (void)(T(2) * m);

    (void)s.to_string();
    (void)v.to_string();
    (void)b.to_string();
    (void)r.to_string();
    (void)m.to_string();

    (void)Multivector<T>(r);
    (void)Rotor<T>(m);
    (void)ga::Pseudoscalar<T>(T(1));
    (void)ga::Complex<T>(s, b);

    // Products: all 25 operator* pairs, plus the named spellings.
    (void)(v * v); (void)(v * b); (void)(v * m); (void)(v * r);
    (void)(b * v); (void)(b * m); (void)(b * r);
    (void)(m * v); (void)(m * b); (void)(m * m); (void)(m * r);
    (void)(r * v); (void)(r * b); (void)(r * m); (void)(r * r);
    (void)(v | v); (void)(v ^ v);
    (void)ga::dot_product(v, v); (void)ga::wedge_product(v, v);
    (void)ga::geometric_product(v, v); (void)ga::geometric_product(m, m);
    (void)ga::rotor_product(v, v); (void)ga::rotor_product(r, r);

    // Sums.
    (void)(s + v); (void)(v + s); (void)(s + b); (void)(b + s); (void)(v + b); (void)(b + v);
    (void)(m + s); (void)(s + m); (void)(m + v); (void)(v + m); (void)(m + b); (void)(b + m);
    (void)(s - v); (void)(v - s); (void)(s - b); (void)(b - s); (void)(v - b); (void)(b - v);
    (void)(m - s); (void)(s - m); (void)(m - v); (void)(v - m); (void)(m - b); (void)(b - m);
    (void)(m + r);
    (void)ga::rotor_sum(s, b); (void)ga::to_rotor(m); (void)ga::to_multivector(r);

    // Grades and contractions.
    (void)ga::grade0(m); (void)ga::grade1(m); (void)ga::grade2(m);
    (void)ga::grade0(r); (void)ga::grade2(r); (void)ga::grade0(s);
    (void)ga::left_contraction(v, v); (void)ga::left_contraction(v, b); (void)ga::left_contraction(b, b);
    (void)ga::right_contraction(v, v); (void)ga::right_contraction(b, v); (void)ga::right_contraction(b, b);
    (void)(v | b); (void)(b | b);
    (void)ga::scalar_product(v, v); (void)ga::scalar_product(b, b);
    (void)ga::scalar_product(r, r); (void)ga::scalar_product(m, m); (void)ga::scalar_product(s, s);

    // Magnitudes, involutions, inverses, duality.
    (void)ga::norm(s); (void)ga::norm(v); (void)ga::norm(b); (void)ga::norm(m); (void)ga::norm(r);
    (void)ga::squared_norm(s); (void)ga::squared_norm(v); (void)ga::squared_norm(b);
    (void)ga::squared_norm(m); (void)ga::squared_norm(r);
    (void)ga::normalize(s); (void)ga::normalize(v); (void)ga::normalize(b);
    (void)ga::normalize(m); (void)ga::normalize(r);
    (void)ga::reverse(s); (void)ga::reverse(v); (void)ga::reverse(b); (void)ga::reverse(m); (void)ga::reverse(r);
    (void)ga::involute(s); (void)ga::involute(v); (void)ga::involute(b);
    (void)ga::involute(m); (void)ga::involute(r);
    (void)ga::conjugate(s); (void)ga::conjugate(v); (void)ga::conjugate(b);
    (void)ga::conjugate(m); (void)ga::conjugate(r);
    (void)ga::inverse(s); (void)ga::inverse(v); (void)ga::inverse(b); (void)ga::inverse(m); (void)ga::inverse(r);
    (void)ga::dual(s); (void)ga::dual(v); (void)ga::dual(b); (void)ga::dual(m);

    // Comparison and printing.
    (void)(s == s); (void)(v == v); (void)(b == b); (void)(m == m); (void)(r == r);
    (void)(s != s); (void)(v != v); (void)(b != b); (void)(m != m); (void)(r != r);
    (void)ga::approx_equal(s, s); (void)ga::approx_equal(v, v); (void)ga::approx_equal(b, b);
    (void)ga::approx_equal(m, m); (void)ga::approx_equal(r, r);
    std::ostringstream os;
    os << s << v << b << m << r;

    // Rotations.
    (void)ga::exp(b); (void)ga::log(r);
    (void)ga::sandwich(v, r); (void)ga::sandwich(b, r); (void)ga::sandwich(m, r);
    (void)ga::rotate(v, r); (void)ga::rotate(b, r); (void)ga::rotate(m, r);
    (void)ga::identity_rotor<T>(); (void)ga::rotor_from_angle(T(1)); (void)ga::rotor_angle(r);
    (void)ga::rotor_between(v, Vector<T>(T(-2), T(1)));
    (void)ga::slerp(r, ga::identity_rotor<T>(), T(0.5));

    // Geometry.
    (void)ga::reflect(v, v); (void)ga::project(v, v); (void)ga::reject(v, v);
    (void)ga::perp(v); (void)ga::signed_area(v, v); (void)ga::angle_between(v, v); (void)ga::angle(v);
}

void test_instantiation_sweep()
{
    section("instantiation sweep");

    instantiate_every_entry_point<float>();
    instantiate_every_entry_point<double>();
    instantiate_every_entry_point<long double>();

    // The sweep is about compiling, not about values; one check records that it
    // ran so a silently deleted call shows up as a missing check.
    check(true, "every entry point instantiated for float, double, long double");
}

} // namespace

int main()
{
    test_construction();
    test_same_type_arithmetic();
    test_bivector_is_the_pseudoscalar();
    test_scalar_products();
    test_compound_assignment();
    test_conversions();
    test_to_string();

    test_basis_multiplication_table();
    test_cayley_table_against_permutation_rule();
    test_multivector_product();
    test_operator_return_types();
    test_mixed_products();
    test_dot_and_wedge();
    test_mixed_grade_addition();
    test_mixed_grade_subtraction();
    test_grade_projection();
    test_contractions();

    test_norms();
    test_normalize();
    test_involutions();
    test_inverse();
    test_dual();
    test_comparison();

    test_exp_log();
    test_rotor_construction();
    test_sandwich();
    test_slerp();
    test_geometry();
    test_complex_isomorphism();

    test_instantiation_sweep();

    std::cout << g_checks << " checks, " << g_failures << " failed.\n";
    if (g_failures != 0) {
        std::cout << "TESTS FAILED\n";
        return 1;
    }
    std::cout << "All CliffordCore Cl2 tests passed\n";
    return 0;
}
