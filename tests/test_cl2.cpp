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

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <type_traits>

#include "../include/cliffordcore/cl2/scalar.hpp"
#include "../include/cliffordcore/cl2/vector.hpp"
#include "../include/cliffordcore/cl2/bivector.hpp"
#include "../include/cliffordcore/cl2/multivector.hpp"
#include "../include/cliffordcore/cl2/rotor.hpp"

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

    test_instantiation_sweep();

    std::cout << g_checks << " checks, " << g_failures << " failed.\n";
    if (g_failures != 0) {
        std::cout << "TESTS FAILED\n";
        return 1;
    }
    std::cout << "All CliffordCore Cl2 tests passed\n";
    return 0;
}
