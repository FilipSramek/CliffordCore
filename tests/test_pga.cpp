// CliffordCore test suite for Cl(3,0,1), namespace CliffordCore::PGA.
//
// Build and run:
//   g++ -std=c++17 -Iinclude tests/test_pga.cpp -o build/test_pga && ./build/test_pga.exe
//   ...or just ./build.sh, which builds every tests/*.cpp
//
// Structure mirrors tests/test_core.cpp: a small assertion harness, then one
// section per type and per operation. Every check records a pass/fail and the
// run continues, so a single build reports every problem rather than stopping
// at the first. This file deliberately includes each header explicitly rather
// than the umbrella, because that is what catches a header that fails to pull
// in its own dependencies.
//
// The sign conventions asserted here are the ones documented in docs/pga.md;
// every one of them was verified against an independent model of the algebra
// before being written down, so a failure means the code, not the expectation.

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <type_traits>

#include "../include/cliffordcore/pga/scalar.hpp"
#include "../include/cliffordcore/pga/vector.hpp"
#include "../include/cliffordcore/pga/bivector.hpp"
#include "../include/cliffordcore/pga/trivector.hpp"
#include "../include/cliffordcore/pga/quadvector.hpp"
#include "../include/cliffordcore/pga/multivector.hpp"
#include "../include/cliffordcore/pga/rotor.hpp"
#include "../include/cliffordcore/pga/translator.hpp"
#include "../include/cliffordcore/pga/motor.hpp"
#include "../include/cliffordcore/pga/operations/dot_product.hpp"
#include "../include/cliffordcore/pga/operations/wedge_product.hpp"
#include "../include/cliffordcore/pga/operations/geometric_product.hpp"
#include "../include/cliffordcore/pga/operations/mixed_products.hpp"
#include "../include/cliffordcore/pga/operations/addition.hpp"
#include "../include/cliffordcore/pga/operations/subtraction.hpp"
#include "../include/cliffordcore/pga/operations/grade.hpp"

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

void check_vector(const Vector<double>& v, double e0, double e1, double e2, double e3,
                  const std::string& what, double tolerance = 1e-12)
{
    check_close(v.e0, e0, what + ".e0", tolerance);
    check_close(v.e1, e1, what + ".e1", tolerance);
    check_close(v.e2, e2, what + ".e2", tolerance);
    check_close(v.e3, e3, what + ".e3", tolerance);
}

void check_bivector(const Bivector<double>& b, double e01, double e02, double e03,
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

void check_trivector(const Trivector<double>& t, double e012, double e013, double e023, double e123,
                     const std::string& what, double tolerance = 1e-12)
{
    check_close(t.e012, e012, what + ".e012", tolerance);
    check_close(t.e013, e013, what + ".e013", tolerance);
    check_close(t.e023, e023, what + ".e023", tolerance);
    check_close(t.e123, e123, what + ".e123", tolerance);
}

void check_quadvector(const Quadvector<double>& q, double e0123, const std::string& what,
                      double tolerance = 1e-12)
{
    check_close(q.e0123, e0123, what + ".e0123", tolerance);
}

void check_rotor(const Rotor<double>& r, double s, double e12, double e13, double e23,
                 const std::string& what, double tolerance = 1e-12)
{
    check_close(r.scalar.value, s, what + ".scalar", tolerance);
    check_close(r.e12, e12, what + ".e12", tolerance);
    check_close(r.e13, e13, what + ".e13", tolerance);
    check_close(r.e23, e23, what + ".e23", tolerance);
}

void check_translator(const Translator<double>& t, double s, double e01, double e02, double e03,
                      const std::string& what, double tolerance = 1e-12)
{
    check_close(t.scalar.value, s, what + ".scalar", tolerance);
    check_close(t.e01, e01, what + ".e01", tolerance);
    check_close(t.e02, e02, what + ".e02", tolerance);
    check_close(t.e03, e03, what + ".e03", tolerance);
}

constexpr double kPi = 3.14159265358979323846;

// Build a multivector from raw components, in grade order and lexicographic
// order within each grade -- the same order as the table in docs/pga.md.
Multivector<double> make_mv(double s,
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
double mv_difference(const Multivector<double>& a, const Multivector<double>& b)
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
void check_mv_close(const Multivector<double>& actual, const Multivector<double>& expected,
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

// ---------------------------------------------------------------------------
// Compile-time checks
//
// Everything in the library is constexpr, so these must evaluate during
// compilation. If constexpr-ness regresses, the build fails rather than the run.
// ---------------------------------------------------------------------------

constexpr Vector<double> kCtA(1.0, 2.0, 3.0, 4.0);
static_assert((kCtA + kCtA).e0 == 2.0, "Vector addition must be constexpr");
static_assert(Scalar<double>(2.0).value == 2.0, "Scalar construction must be constexpr");
static_assert(Vector<double>().e1 == 0.0, "default Vector must be constexpr and zeroed");
static_assert(Bivector<double>(1, 2, 3, 4, 5, 6).e23 == 6.0, "Bivector construction must be constexpr");
static_assert(Quadvector<double>(3.0).e0123 == 3.0, "Quadvector construction must be constexpr");

// Widening is implicit, narrowing is not. These pin the conversion design so
// that `motor + rotor` keeps compiling and `rotor + multivector` keeps failing.
static_assert(std::is_convertible<Rotor<double>, Motor<double>>::value, "Rotor widens to Motor implicitly");
static_assert(std::is_convertible<Translator<double>, Motor<double>>::value, "Translator widens to Motor implicitly");
static_assert(std::is_convertible<Rotor<double>, Multivector<double>>::value, "Rotor widens to Multivector implicitly");
static_assert(std::is_convertible<Translator<double>, Multivector<double>>::value, "Translator widens to Multivector implicitly");
static_assert(std::is_convertible<Motor<double>, Multivector<double>>::value, "Motor widens to Multivector implicitly");
static_assert(!std::is_convertible<Multivector<double>, Motor<double>>::value, "Multivector does not narrow to Motor implicitly");
static_assert(!std::is_convertible<Multivector<double>, Rotor<double>>::value, "Multivector does not narrow to Rotor implicitly");
static_assert(!std::is_convertible<Multivector<double>, Translator<double>>::value, "Multivector does not narrow to Translator implicitly");
static_assert(!std::is_convertible<Motor<double>, Rotor<double>>::value, "Motor does not narrow to Rotor implicitly");
static_assert(!std::is_convertible<Motor<double>, Translator<double>>::value, "Motor does not narrow to Translator implicitly");
static_assert(std::is_constructible<Rotor<double>, Multivector<double>>::value, "Rotor narrows from Multivector explicitly");
static_assert(std::is_constructible<Motor<double>, Multivector<double>>::value, "Motor narrows from Multivector explicitly");

// The aliases are the same types, not new ones.
static_assert(std::is_same<ga::Plane<double>, Vector<double>>::value, "Plane is Vector");
static_assert(std::is_same<ga::Line<double>, Bivector<double>>::value, "Line is Bivector");
static_assert(std::is_same<ga::Twist<double>, Bivector<double>>::value, "Twist is Bivector");
static_assert(std::is_same<ga::Point<double>, Trivector<double>>::value, "Point is Trivector");

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

void test_scalar()
{
    section("Scalar");

    check_scalar(Scalar<double>(), 0.0, "default ctor zeroes");
    check_scalar(Scalar<double>(3.5), 3.5, "value ctor");

    const Scalar<double> a(6.0);
    const Scalar<double> b(4.0);

    check_scalar(a + b, 10.0, "operator+");
    check_scalar(a - b, 2.0, "operator-");
    check_scalar(a * b, 24.0, "operator* (Scalar)");
    check_scalar(a / b, 1.5, "operator/");
    check_scalar(-a, -6.0, "unary operator-");

    // Cross-type multiplication: Scalar on the left, every other type on the right.
    const Scalar<double> two(2.0);
    check_vector(two * Vector<double>(1.0, 2.0, 3.0, 4.0), 2.0, 4.0, 6.0, 8.0, "Scalar * Vector");
    check_bivector(two * Bivector<double>(1.0, 2.0, 3.0, 4.0, 5.0, 6.0),
                   2.0, 4.0, 6.0, 8.0, 10.0, 12.0, "Scalar * Bivector");
    check_trivector(two * Trivector<double>(1.0, 2.0, 3.0, 4.0), 2.0, 4.0, 6.0, 8.0, "Scalar * Trivector");
    check_quadvector(two * Quadvector<double>(5.0), 10.0, "Scalar * Quadvector");

    const Multivector<double> mv = make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
    check_mv_close(two * mv, make_mv(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32),
                   "Scalar * Multivector distributes over every grade");

    check_rotor(two * Rotor<double>(Scalar<double>(1.0), 2.0, 3.0, 4.0), 2.0, 4.0, 6.0, 8.0, "Scalar * Rotor");
    check_translator(two * Translator<double>(Scalar<double>(1.0), 2.0, 3.0, 4.0), 2.0, 4.0, 6.0, 8.0,
                     "Scalar * Translator");

    const Motor<double> mo(Scalar<double>(1.0), Bivector<double>(1, 2, 3, 4, 5, 6), Quadvector<double>(7.0));
    const Motor<double> scaledMo = two * mo;
    check_scalar(scaledMo.scalar, 2.0, "Scalar * Motor scalar");
    check_bivector(scaledMo.bivector, 2, 4, 6, 8, 10, 12, "Scalar * Motor bivector");
    check_quadvector(scaledMo.quadvector, 14.0, "Scalar * Motor quadvector");
}

void test_vector()
{
    section("Vector");

    check_vector(Vector<double>(), 0, 0, 0, 0, "default ctor zeroes");
    check_vector(Vector<double>(1, 2, 3, 4), 1, 2, 3, 4, "component ctor");

    const Vector<double> a(1, 2, 3, 4);
    const Vector<double> b(5, 6, 7, 8);

    check_vector(a + b, 6, 8, 10, 12, "operator+");
    check_vector(a - b, -4, -4, -4, -4, "operator-");
    check_vector(-a, -1, -2, -3, -4, "unary operator-");
    check_vector(a * Scalar<double>(2.0), 2, 4, 6, 8, "operator* (Scalar)");
    check_vector(a / Scalar<double>(2.0), 0.5, 1, 1.5, 2, "operator/ (Scalar)");
    check_vector(a * 2.0, 2, 4, 6, 8, "operator* (raw T via implicit Scalar)");
    check_vector(2.0 * a, 2, 4, 6, 8, "raw T on the left");

    // The magnitude is Euclidean: e0 does not contribute, because e0^2 = 0.
    check_scalar(a.magnitude(), std::sqrt(4.0 + 9.0 + 16.0), "magnitude ignores e0");
    check_scalar(Vector<double>(7, 0, 0, 0).magnitude(), 0.0, "an ideal plane has zero magnitude");
}

void test_bivector()
{
    section("Bivector");

    check_bivector(Bivector<double>(), 0, 0, 0, 0, 0, 0, "default ctor zeroes");
    check_bivector(Bivector<double>(1, 2, 3, 4, 5, 6), 1, 2, 3, 4, 5, 6, "component ctor");

    const Bivector<double> a(1, 2, 3, 4, 5, 6);
    const Bivector<double> b(6, 5, 4, 3, 2, 1);

    check_bivector(a + b, 7, 7, 7, 7, 7, 7, "operator+");
    check_bivector(a - b, -5, -3, -1, 1, 3, 5, "operator-");
    check_bivector(-a, -1, -2, -3, -4, -5, -6, "unary operator-");
    check_bivector(a * Scalar<double>(2.0), 2, 4, 6, 8, 10, 12, "operator* (Scalar)");
    check_bivector(a / Scalar<double>(2.0), 0.5, 1, 1.5, 2, 2.5, 3, "operator/ (Scalar)");
    check_bivector(2.0 * a, 2, 4, 6, 8, 10, 12, "raw T on the left");

    // Euclidean magnitude: the direction part only.
    check_scalar(a.magnitude(), std::sqrt(16.0 + 25.0 + 36.0), "magnitude ignores the ideal part");
    check_scalar(Bivector<double>(1, 2, 3, 0, 0, 0).magnitude(), 0.0, "an ideal line has zero magnitude");
}

void test_trivector()
{
    section("Trivector");

    check_trivector(Trivector<double>(), 0, 0, 0, 0, "default ctor zeroes");
    check_trivector(Trivector<double>(1, 2, 3, 4), 1, 2, 3, 4, "component ctor");

    const Trivector<double> a(1, 2, 3, 4);
    const Trivector<double> b(4, 3, 2, 1);

    check_trivector(a + b, 5, 5, 5, 5, "operator+");
    check_trivector(a - b, -3, -1, 1, 3, "operator-");
    check_trivector(-a, -1, -2, -3, -4, "unary operator-");
    check_trivector(a * Scalar<double>(2.0), 2, 4, 6, 8, "operator* (Scalar)");
    check_trivector(a / Scalar<double>(2.0), 0.5, 1, 1.5, 2, "operator/ (Scalar)");
    check_trivector(2.0 * a, 2, 4, 6, 8, "raw T on the left");

    // The weight is the only Euclidean component.
    check_scalar(a.magnitude(), 4.0, "magnitude is the absolute weight");
    check_scalar(Trivector<double>(1, 2, 3, -4).magnitude(), 4.0, "magnitude is non-negative");
    check_scalar(Trivector<double>(1, 2, 3, 0).magnitude(), 0.0, "an ideal point has zero magnitude");
}

void test_quadvector()
{
    section("Quadvector");

    check_quadvector(Quadvector<double>(), 0.0, "default ctor zeroes");
    check_quadvector(Quadvector<double>(7.0), 7.0, "value ctor");

    const Quadvector<double> a(6.0);
    const Quadvector<double> b(3.0);

    check_quadvector(a + b, 9.0, "operator+");
    check_quadvector(a - b, 3.0, "operator-");
    check_quadvector(-a, -6.0, "unary operator-");

    // Two pseudoscalars multiply to zero, always: e0123 contains e0 twice
    // once squared. The return type is Scalar so that every pair has an
    // operator*, but the value is identically zero.
    static_assert(std::is_same<decltype(a * b), Scalar<double>>::value,
                  "Quadvector * Quadvector produces a Scalar");
    check((a * b).value == 0.0, "e0123 * e0123 == 0 exactly");

    check_quadvector(a * Scalar<double>(2.0), 12.0, "operator* (Scalar)");
    check_quadvector(a / Scalar<double>(2.0), 3.0, "operator/ (Scalar)");

    // Raw numeric overloads exist so `q * 2.0` stays unambiguous now that both
    // Quadvector and Scalar are constructible from one T.
    check_quadvector(a * 2.0, 12.0, "operator* (raw T)");
    check_quadvector(a / 2.0, 3.0, "operator/ (raw T)");
    check_quadvector(2.0 * a, 12.0, "raw T on the left");
}

void test_multivector()
{
    section("Multivector");

    const Multivector<double> zero;
    check_mv_close(zero, make_mv(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), "default ctor zeroes");

    const Multivector<double> a = make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
    const Multivector<double> b = make_mv(16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1);

    check_mv_close(a + b, make_mv(17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17), "operator+");
    check_mv_close(a - b, make_mv(-15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15), "operator-");
    check_mv_close(-a, make_mv(-1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15, -16), "unary operator-");
    check_mv_close(a * Scalar<double>(2.0), make_mv(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32),
                   "operator* (Scalar)");
    check_mv_close(a / Scalar<double>(2.0), make_mv(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8),
                   "operator/ (Scalar)");
    check_mv_close(2.0 * a, a * Scalar<double>(2.0), "raw T on the left");
    check_mv_close(zero - a, -a, "0 - a == -a");

    // Widening from the even types is implicit and lossless.
    const Rotor<double> r(Scalar<double>(1.0), 2.0, 3.0, 4.0);
    const Multivector<double> fromRotor = r;
    check_mv_close(fromRotor, make_mv(1, 0, 0, 0, 0, 0, 0, 0, 2, 3, 4, 0, 0, 0, 0, 0),
                   "Multivector(Rotor) fills the Euclidean bivector slots");

    const Translator<double> t(Scalar<double>(1.0), 2.0, 3.0, 4.0);
    const Multivector<double> fromTranslator = t;
    check_mv_close(fromTranslator, make_mv(1, 0, 0, 0, 0, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0),
                   "Multivector(Translator) fills the ideal bivector slots");

    const Motor<double> mo(Scalar<double>(1.0), Bivector<double>(2, 3, 4, 5, 6, 7), Quadvector<double>(8.0));
    const Multivector<double> fromMotor = mo;
    check_mv_close(fromMotor, make_mv(1, 0, 0, 0, 0, 2, 3, 4, 5, 6, 7, 0, 0, 0, 0, 8),
                   "Multivector(Motor) fills grades 0, 2 and 4");
}

void test_rotor()
{
    section("Rotor");

    check_rotor(Rotor<double>(), 0, 0, 0, 0, "default ctor zeroes");
    check_rotor(Rotor<double>(Scalar<double>(1.0), 2.0, 3.0, 4.0), 1, 2, 3, 4, "component ctor");

    const Rotor<double> a(Scalar<double>(1.0), 2.0, 3.0, 4.0);
    const Rotor<double> b(Scalar<double>(4.0), 3.0, 2.0, 1.0);

    check_rotor(a + b, 5, 5, 5, 5, "operator+");
    check_rotor(a - b, -3, -1, 1, 3, "operator-");
    check_rotor(-a, -1, -2, -3, -4, "unary operator-");
    check_rotor(a * Scalar<double>(2.0), 2, 4, 6, 8, "operator* (Scalar)");
    check_rotor(a / Scalar<double>(2.0), 0.5, 1, 1.5, 2, "operator/ (Scalar)");
    check_rotor(2.0 * a, 2, 4, 6, 8, "raw T on the left");

    // Narrowing is explicit and keeps only the rotor slots.
    const Multivector<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
    check_rotor(Rotor<double>(m), 1, 9, 10, 11, "explicit Rotor(Multivector) keeps scalar and Euclidean bivector");

    const Motor<double> mo(Scalar<double>(1.0), Bivector<double>(2, 3, 4, 5, 6, 7), Quadvector<double>(8.0));
    check_rotor(Rotor<double>(mo), 1, 5, 6, 7, "explicit Rotor(Motor) drops the ideal part and e0123");
}

void test_translator()
{
    section("Translator");

    check_translator(Translator<double>(), 0, 0, 0, 0, "default ctor zeroes");
    check_translator(Translator<double>(Scalar<double>(1.0), 2.0, 3.0, 4.0), 1, 2, 3, 4, "component ctor");

    const Translator<double> a(Scalar<double>(1.0), 2.0, 3.0, 4.0);
    const Translator<double> b(Scalar<double>(4.0), 3.0, 2.0, 1.0);

    check_translator(a + b, 5, 5, 5, 5, "operator+");
    check_translator(a - b, -3, -1, 1, 3, "operator-");
    check_translator(-a, -1, -2, -3, -4, "unary operator-");
    check_translator(a * Scalar<double>(2.0), 2, 4, 6, 8, "operator* (Scalar)");
    check_translator(a / Scalar<double>(2.0), 0.5, 1, 1.5, 2, "operator/ (Scalar)");
    check_translator(2.0 * a, 2, 4, 6, 8, "raw T on the left");

    const Multivector<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
    check_translator(Translator<double>(m), 1, 6, 7, 8, "explicit Translator(Multivector) keeps scalar and ideal bivector");

    const Motor<double> mo(Scalar<double>(1.0), Bivector<double>(2, 3, 4, 5, 6, 7), Quadvector<double>(8.0));
    check_translator(Translator<double>(mo), 1, 2, 3, 4, "explicit Translator(Motor) drops the Euclidean part and e0123");
}

void test_motor()
{
    section("Motor");

    const Motor<double> zero;
    check_scalar(zero.scalar, 0.0, "default ctor scalar");
    check_bivector(zero.bivector, 0, 0, 0, 0, 0, 0, "default ctor bivector");
    check_quadvector(zero.quadvector, 0.0, "default ctor quadvector");

    const Motor<double> a(Scalar<double>(1.0), Bivector<double>(2, 3, 4, 5, 6, 7), Quadvector<double>(8.0));
    const Motor<double> b(Scalar<double>(8.0), Bivector<double>(7, 6, 5, 4, 3, 2), Quadvector<double>(1.0));

    const Motor<double> sum = a + b;
    check_scalar(sum.scalar, 9.0, "operator+ scalar");
    check_bivector(sum.bivector, 9, 9, 9, 9, 9, 9, "operator+ bivector");
    check_quadvector(sum.quadvector, 9.0, "operator+ quadvector");

    const Motor<double> diff = a - b;
    check_scalar(diff.scalar, -7.0, "operator- scalar");
    check_bivector(diff.bivector, -5, -3, -1, 1, 3, 5, "operator- bivector");
    check_quadvector(diff.quadvector, 7.0, "operator- quadvector");

    const Motor<double> neg = -a;
    check_scalar(neg.scalar, -1.0, "unary operator- scalar");
    check_bivector(neg.bivector, -2, -3, -4, -5, -6, -7, "unary operator- bivector");
    check_quadvector(neg.quadvector, -8.0, "unary operator- quadvector");

    const Motor<double> scaled = a * Scalar<double>(2.0);
    check_scalar(scaled.scalar, 2.0, "operator* (Scalar) scalar");
    check_bivector(scaled.bivector, 4, 6, 8, 10, 12, 14, "operator* (Scalar) bivector");
    check_quadvector(scaled.quadvector, 16.0, "operator* (Scalar) quadvector");

    const Motor<double> halved = a / Scalar<double>(2.0);
    check_scalar(halved.scalar, 0.5, "operator/ (Scalar) scalar");
    check_bivector(halved.bivector, 1, 1.5, 2, 2.5, 3, 3.5, "operator/ (Scalar) bivector");
    check_quadvector(halved.quadvector, 4.0, "operator/ (Scalar) quadvector");

    const Motor<double> left = 2.0 * a;
    check_scalar(left.scalar, 2.0, "raw T on the left scalar");
    check_quadvector(left.quadvector, 16.0, "raw T on the left quadvector");

    // Widening from Rotor and Translator is implicit and lands in the right slots.
    const Motor<double> fromRotor = Rotor<double>(Scalar<double>(1.0), 2.0, 3.0, 4.0);
    check_scalar(fromRotor.scalar, 1.0, "Motor(Rotor) scalar");
    check_bivector(fromRotor.bivector, 0, 0, 0, 2, 3, 4, "Motor(Rotor) fills the Euclidean bivector slots");
    check_quadvector(fromRotor.quadvector, 0.0, "Motor(Rotor) has no e0123");

    const Motor<double> fromTranslator = Translator<double>(Scalar<double>(1.0), 2.0, 3.0, 4.0);
    check_scalar(fromTranslator.scalar, 1.0, "Motor(Translator) scalar");
    check_bivector(fromTranslator.bivector, 2, 3, 4, 0, 0, 0, "Motor(Translator) fills the ideal bivector slots");
    check_quadvector(fromTranslator.quadvector, 0.0, "Motor(Translator) has no e0123");

    // Narrowing from Multivector is explicit and keeps the even grades only.
    const Multivector<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
    const Motor<double> fromMv(m);
    check_scalar(fromMv.scalar, 1.0, "explicit Motor(Multivector) scalar");
    check_bivector(fromMv.bivector, 6, 7, 8, 9, 10, 11, "explicit Motor(Multivector) bivector");
    check_quadvector(fromMv.quadvector, 16.0, "explicit Motor(Multivector) quadvector");

    // A motor plus a rotor compiles through the implicit widening; the reverse
    // order does not (a member operator's left operand never converts), which
    // is the same asymmetry Cl3 documents for multivector + rotor.
    const Motor<double> mixed = a + Rotor<double>(Scalar<double>(1.0), 1.0, 1.0, 1.0);
    check_bivector(mixed.bivector, 2, 3, 4, 6, 7, 8, "motor + rotor widens the rotor");
}

void test_compound_assignment()
{
    section("compound assignment");

    { Scalar<double> s(2.0); s += Scalar<double>(1.0); s -= Scalar<double>(0.5); s *= Scalar<double>(4.0); s /= Scalar<double>(2.0);
      check_scalar(s, 5.0, "Scalar += -= *= /="); }
    { Vector<double> v(1, 2, 3, 4); v += Vector<double>(1, 1, 1, 1); v -= Vector<double>(0, 1, 0, 1); v *= 2.0; v /= Scalar<double>(2.0);
      check_vector(v, 2, 2, 4, 4, "Vector += -= *= /="); }
    { Bivector<double> b(1, 2, 3, 4, 5, 6); b += b; b -= Bivector<double>(1, 1, 1, 1, 1, 1); b *= Scalar<double>(2.0); b /= 2.0;
      check_bivector(b, 1, 3, 5, 7, 9, 11, "Bivector += -= *= /="); }
    { Trivector<double> t(1, 2, 3, 4); t += t; t -= Trivector<double>(1, 1, 1, 1); t *= 2.0; t /= Scalar<double>(2.0);
      check_trivector(t, 1, 3, 5, 7, "Trivector += -= *= /="); }
    { Quadvector<double> q(3.0); q += q; q -= Quadvector<double>(1.0); q *= 2.0; q /= Scalar<double>(2.0);
      check_quadvector(q, 5.0, "Quadvector += -= *= /="); }
    { Multivector<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
      m += m; m -= make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16); m *= Scalar<double>(2.0); m /= Scalar<double>(2.0);
      check_mv_close(m, make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), "Multivector += -= *= /="); }
    { Rotor<double> r(Scalar<double>(1.0), 2.0, 3.0, 4.0); r += r; r -= Rotor<double>(Scalar<double>(1.0), 1, 1, 1); r *= Scalar<double>(2.0); r /= Scalar<double>(2.0);
      check_rotor(r, 1, 3, 5, 7, "Rotor += -= *= /="); }
    { Translator<double> t(Scalar<double>(1.0), 2.0, 3.0, 4.0); t += t; t -= Translator<double>(Scalar<double>(1.0), 1, 1, 1); t *= Scalar<double>(2.0); t /= Scalar<double>(2.0);
      check_translator(t, 1, 3, 5, 7, "Translator += -= *= /="); }
    { Motor<double> m(Scalar<double>(1.0), Bivector<double>(1, 1, 1, 1, 1, 1), Quadvector<double>(1.0));
      m += m; m -= Motor<double>(Scalar<double>(1.0), Bivector<double>(), Quadvector<double>()); m *= Scalar<double>(2.0); m /= Scalar<double>(2.0);
      check_scalar(m.scalar, 1.0, "Motor += -= *= /= scalar");
      check_bivector(m.bivector, 2, 2, 2, 2, 2, 2, "Motor += -= *= /= bivector");
      check_quadvector(m.quadvector, 2.0, "Motor += -= *= /= quadvector"); }
}

// ---------------------------------------------------------------------------
// Products
// ---------------------------------------------------------------------------

// The sixteen basis blades in lexicographic order, as index lists. Index 0 is
// the scalar (empty list). This order is the one make_mv() takes.
struct BladeIndices {
    int count;
    int index[4];
};

const BladeIndices kBlades[16] = {
    {0, {0, 0, 0, 0}},
    {1, {0, 0, 0, 0}}, {1, {1, 0, 0, 0}}, {1, {2, 0, 0, 0}}, {1, {3, 0, 0, 0}},
    {2, {0, 1, 0, 0}}, {2, {0, 2, 0, 0}}, {2, {0, 3, 0, 0}},
    {2, {1, 2, 0, 0}}, {2, {1, 3, 0, 0}}, {2, {2, 3, 0, 0}},
    {3, {0, 1, 2, 0}}, {3, {0, 1, 3, 0}}, {3, {0, 2, 3, 0}}, {3, {1, 2, 3, 0}},
    {4, {0, 1, 2, 3}},
};

const char* const kBladeNames[16] = {
    "1", "e0", "e1", "e2", "e3", "e01", "e02", "e03", "e12", "e13", "e23",
    "e012", "e013", "e023", "e123", "e0123",
};

// A multivector with a single basis blade set to one.
Multivector<double> basis_blade(int which)
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
int reference_blade_product(int i, int j, int& result_blade)
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

// Local stand-ins for the primitives that arrive with primitives.hpp; written
// straight from docs/pga.md so the products can be tested on real geometry.
Trivector<double> test_point(double x, double y, double z)
{
    return Trivector<double>(-z, y, -x, 1.0);
}

Vector<double> test_plane(double a, double b, double c, double d)
{
    return Vector<double>(d, a, b, c);
}

void test_basis_multiplication_table()
{
    section("basis table");

    const Multivector<double> e0 = basis_blade(1), e1 = basis_blade(2), e2 = basis_blade(3), e3 = basis_blade(4);
    const Multivector<double> zero;

    // The metric: e0 is null, the rest square to +1.
    check(mv_difference(e0 * e0, zero) == 0.0, "e0 * e0 == 0 exactly");
    check_mv_close(e1 * e1, basis_blade(0), "e1 * e1 == 1");
    check_mv_close(e2 * e2, basis_blade(0), "e2 * e2 == 1");
    check_mv_close(e3 * e3, basis_blade(0), "e3 * e3 == 1");

    // Anticommutation of distinct basis vectors, e0 included.
    check_mv_close(e0 * e1, -(e1 * e0), "e0 and e1 anticommute");
    check_mv_close(e1 * e2, -(e2 * e1), "e1 and e2 anticommute");

    // The spot checks from docs/pga.md.
    check_mv_close(basis_blade(5) * basis_blade(8), basis_blade(6), "e01 * e12 == +e02");
    check_mv_close(basis_blade(8) * basis_blade(5), -basis_blade(6), "e12 * e01 == -e02");
    check_mv_close(basis_blade(14) * basis_blade(15), basis_blade(1), "e123 * e0123 == +e0");
    check_mv_close(basis_blade(15) * basis_blade(14), -basis_blade(1), "e0123 * e123 == -e0");
    check_mv_close(basis_blade(11) * basis_blade(14), -basis_blade(7), "e012 * e123 == -e03");
    check(mv_difference(basis_blade(11) * basis_blade(12), zero) == 0.0, "e012 * e013 == 0 exactly");
    check(mv_difference(basis_blade(15) * basis_blade(15), zero) == 0.0, "e0123 * e0123 == 0 exactly");
    check_mv_close(basis_blade(14) * basis_blade(14), -basis_blade(0), "e123 * e123 == -1");
    check_mv_close(basis_blade(8) * basis_blade(15), -basis_blade(7), "e12 * e0123 == -e03");
    check_mv_close(basis_blade(9) * basis_blade(15), basis_blade(6), "e13 * e0123 == +e02");
    check_mv_close(basis_blade(10) * basis_blade(15), -basis_blade(5), "e23 * e0123 == -e01");

    // The pseudoscalar commutes with every bivector and anticommutes with every vector.
    const Multivector<double> I = basis_blade(15);
    for (int b = 5; b <= 10; ++b) {
        check_mv_close(I * basis_blade(b), basis_blade(b) * I, std::string("e0123 commutes with ") + kBladeNames[b]);
    }
    for (int v = 1; v <= 4; ++v) {
        check_mv_close(I * basis_blade(v), -(basis_blade(v) * I), std::string("e0123 anticommutes with ") + kBladeNames[v]);
    }
}

void test_cayley_table_against_permutation_rule()
{
    section("Cayley table vs permutation rule");

    // Every one of the 256 basis products, against the independent reference.
    // One check per entry keeps a failure pointing at the exact term.
    int nonzero = 0;
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            int k = -1;
            const int sign = reference_blade_product(i, j, k);
            const Multivector<double> product = basis_blade(i) * basis_blade(j);
            const std::string what = std::string(kBladeNames[i]) + " * " + kBladeNames[j];
            if (sign == 0) {
                check(mv_difference(product, Multivector<double>()) == 0.0, what + " == 0");
            } else {
                ++nonzero;
                check_mv_close(product, static_cast<double>(sign) * basis_blade(k),
                               what + " == " + (sign > 0 ? "+" : "-") + kBladeNames[k]);
            }
        }
    }
    check(nonzero == 192, "192 of the 256 basis products are non-zero");
}

void test_multivector_product()
{
    section("multivector product");

    const Multivector<double> a = make_mv(1, -2, 3, 0.5, -1, 2, -3, 1, 0.25, -0.5, 4, 1, -1, 2, -2, 0.75);
    const Multivector<double> b = make_mv(-0.5, 1, 2, -3, 4, -1, 0.5, 2, 3, -2, 1, 0.5, 1.5, -1, 1, -2);
    const Multivector<double> c = make_mv(2, 0.5, -1, 1, 1, 3, -2, 0.5, -1, 1, -0.5, -2, 1, 0.5, 3, 1);
    const Multivector<double> one = basis_blade(0);

    check_mv_close((a * b) * c, a * (b * c), "associativity", 1e-9);
    check_mv_close(a * (b + c), a * b + a * c, "left distributivity", 1e-9);
    check_mv_close((a + b) * c, a * c + b * c, "right distributivity", 1e-9);
    check_mv_close(one * a, a, "1 * a == a");
    check_mv_close(a * one, a, "a * 1 == a");
    check(mv_difference(a * b, b * a) > 1e-6, "the product is not commutative");

    // operator* and the named function are the same thing.
    check_mv_close(a * b, ga::geometric_product(a, b), "operator* == geometric_product");
}

void test_dot_and_wedge()
{
    section("dot and wedge");

    const Vector<double> a(1, 2, 3, 4);   // e0, e1, e2, e3
    const Vector<double> b(5, 6, 7, 8);

    // The dot product ignores e0: only the normals matter.
    check_scalar(a | b, 2 * 6 + 3 * 7 + 4 * 8, "operator| ignores e0");
    check_scalar(ga::dot_product(a, b), (a | b).value, "dot_product == operator|");
    check_scalar(Vector<double>(1, 0, 0, 0) | Vector<double>(1, 0, 0, 0), 0.0, "e0 | e0 == 0");
    check_scalar(a | b, (b | a).value, "dot product is symmetric");

    // The wedge is antisymmetric and lexicographic.
    check_bivector(a ^ b, 1 * 6 - 2 * 5, 1 * 7 - 3 * 5, 1 * 8 - 4 * 5,
                   2 * 7 - 3 * 6, 2 * 8 - 4 * 6, 3 * 8 - 4 * 7, "operator^ components");
    check_bivector(b ^ a, -(a ^ b).e01, -(a ^ b).e02, -(a ^ b).e03, -(a ^ b).e12, -(a ^ b).e13, -(a ^ b).e23,
                   "wedge is antisymmetric");
    const Bivector<double> self = a ^ a;
    check(self.e01 == 0.0 && self.e02 == 0.0 && self.e03 == 0.0
          && self.e12 == 0.0 && self.e13 == 0.0 && self.e23 == 0.0, "v ^ v == 0 exactly");
    check_bivector(ga::wedge_product(a, b), (a ^ b).e01, (a ^ b).e02, (a ^ b).e03, (a ^ b).e12, (a ^ b).e13, (a ^ b).e23,
                   "wedge_product == operator^");

    // Both agree with the general product.
    const Multivector<double> product = a * b;
    check_scalar(product.scalar, (a | b).value, "scalar part of v*w is the dot");
    check_bivector(product.bivector, (a ^ b).e01, (a ^ b).e02, (a ^ b).e03, (a ^ b).e12, (a ^ b).e13, (a ^ b).e23,
                   "bivector part of v*w is the wedge");
    check_mv_close(ga::geometric_product(a, b), ga::geometric_product(ga::detail::promote(a), ga::detail::promote(b)),
                   "Vector*Vector agrees with the 16x16 table");

    // Geometry: the wedge of two planes is the line where they meet.
    const Vector<double> x0 = test_plane(1, 0, 0, 0);   // x = 0
    const Vector<double> y0 = test_plane(0, 1, 0, 0);   // y = 0
    const Vector<double> x1 = test_plane(1, 0, 0, -1);  // x = 1
    check_bivector(x0 ^ y0, 0, 0, 0, 1, 0, 0, "meet(x=0, y=0) is the z-axis e12");
    check_bivector(x1 ^ y0, 0, -1, 0, 1, 0, 0, "meet(x=1, y=0) is e12 - e02");
    const Vector<double> z0 = test_plane(0, 0, 1, 0);
    const Vector<double> z1 = test_plane(0, 0, 1, -1);
    check_scalar((z0 ^ z1).magnitude(), 0.0, "parallel planes meet in an ideal line");
    check_bivector(z0 ^ z1, 0, 0, 1, 0, 0, 0, "meet(z=0, z=1) is the ideal line +e03");
}

void test_incidence()
{
    section("incidence");

    // plane ^ point = (a x + b y + c z + d) e0123, verified in docs/pga.md.
    const Vector<double> p = test_plane(1, 2, 3, 4);
    const Trivector<double> P = test_point(0.5, -1, 2);
    check_quadvector(p ^ P, 1 * 0.5 + 2 * -1 + 3 * 2 + 4, "plane ^ point is the plane equation");
    check_quadvector(P ^ p, -(1 * 0.5 + 2 * -1 + 3 * 2 + 4), "point ^ plane has the opposite sign");

    // A point on the plane gives exactly zero.
    const Vector<double> x1 = test_plane(1, 0, 0, -1);
    check((x1 ^ test_point(1, 7, -3)).e0123 == 0.0, "point(1, 7, -3) lies on x = 1");
    check((x1 ^ test_point(2, 0, 0)).e0123 != 0.0, "point(2, 0, 0) does not lie on x = 1");

    // The origin is e123 and lies on every plane through the origin.
    check((test_plane(3, -2, 5, 0) ^ test_point(0, 0, 0)).e0123 == 0.0, "the origin lies on 3x - 2y + 5z = 0");

    // The point is the triple meet of its coordinate planes.
    const Vector<double> px = test_plane(1, 0, 0, -1), py = test_plane(0, 1, 0, -2), pz = test_plane(0, 0, 1, -3);
    check_trivector((px ^ py) ^ pz, -3, 2, -1, 1, "point(1,2,3) == (x=1) ^ (y=2) ^ (z=3)");
    check_trivector(px ^ (py ^ pz), -3, 2, -1, 1, "the triple meet associates");

    // Two lines: zero iff coplanar.
    const Bivector<double> zAxis(0, 0, 0, 1, 0, 0);
    const Bivector<double> zParallel(0, -1, 0, 1, 0, 0);   // through (1,0,0)
    const Bivector<double> xAxis(0, 0, 0, 0, 0, 1);
    check((zAxis ^ zParallel).e0123 == 0.0, "parallel lines are coplanar");
    check((zAxis ^ xAxis).e0123 == 0.0, "intersecting lines are coplanar");
    const Bivector<double> skew = test_plane(0, 1, 0, -1) ^ test_plane(0, 0, 1, 0);   // y = 1, z = 0: parallel to x, offset
    check((zAxis ^ skew).e0123 != 0.0, "skew lines have a non-zero wedge");
}

void test_even_closure()
{
    section("even subalgebra closure");

    const Rotor<double> r(Scalar<double>(0.5), 0.25, -0.75, 1.0);
    const Rotor<double> s(Scalar<double>(-1.0), 2.0, 0.5, -0.25);
    const Translator<double> t(Scalar<double>(1.0), 0.5, -1.0, 2.0);
    const Translator<double> u(Scalar<double>(2.0), -0.25, 3.0, 1.0);
    const Motor<double> m(Scalar<double>(0.75), Bivector<double>(1, -2, 0.5, 3, -1, 2), Quadvector<double>(-0.5));
    const Motor<double> n(Scalar<double>(-1.5), Bivector<double>(0.5, 1, -1, -2, 0.25, 1), Quadvector<double>(2.0));

    // The closed forms must agree with the one multiplication table.
    check_mv_close(ga::to_multivector(r * s), ga::to_multivector(r) * ga::to_multivector(s), "Rotor * Rotor matches the table");
    check_mv_close(ga::to_multivector(t * u), ga::to_multivector(t) * ga::to_multivector(u), "Translator * Translator matches the table");
    check_mv_close(ga::to_multivector(m * n), ga::to_multivector(m) * ga::to_multivector(n), "Motor * Motor matches the table");
    check_mv_close(ga::to_multivector(r * t), ga::to_multivector(r) * ga::to_multivector(t), "Rotor * Translator matches the table");
    check_mv_close(ga::to_multivector(t * r), ga::to_multivector(t) * ga::to_multivector(r), "Translator * Rotor matches the table");
    check_mv_close(ga::to_multivector(m * r), ga::to_multivector(m) * ga::to_multivector(r), "Motor * Rotor matches the table");
    check_mv_close(ga::to_multivector(r * m), ga::to_multivector(r) * ga::to_multivector(m), "Rotor * Motor matches the table");
    check_mv_close(ga::to_multivector(m * t), ga::to_multivector(m) * ga::to_multivector(t), "Motor * Translator matches the table");
    check_mv_close(ga::to_multivector(t * m), ga::to_multivector(t) * ga::to_multivector(m), "Translator * Motor matches the table");

    // Return types: closures stay narrow, everything else is a Multivector.
    static_assert(std::is_same<decltype(r * s), Rotor<double>>::value, "Rotor * Rotor is a Rotor");
    static_assert(std::is_same<decltype(t * u), Translator<double>>::value, "Translator * Translator is a Translator");
    static_assert(std::is_same<decltype(m * n), Motor<double>>::value, "Motor * Motor is a Motor");
    static_assert(std::is_same<decltype(r * t), Motor<double>>::value, "Rotor * Translator is a Motor");
    static_assert(std::is_same<decltype(t * r), Motor<double>>::value, "Translator * Rotor is a Motor");
    static_assert(std::is_same<decltype(m * r), Motor<double>>::value, "Motor * Rotor is a Motor");
    static_assert(std::is_same<decltype(t * m), Motor<double>>::value, "Translator * Motor is a Motor");

    // The named spellings are the same products.
    check_mv_close(ga::to_multivector(ga::rotor_product(r, s)), ga::to_multivector(r * s), "rotor_product == operator*");
    check_mv_close(ga::to_multivector(ga::translator_product(t, u)), ga::to_multivector(t * u), "translator_product == operator*");
    check_mv_close(ga::to_multivector(ga::motor_product(m, n)), ga::to_multivector(m * n), "motor_product == operator*");

    // Normalised translators compose by adding their displacements.
    const Translator<double> tx(Scalar<double>(1.0), -0.5, 0, 0);   // +1 along x
    const Translator<double> ty(Scalar<double>(1.0), 0, -1.0, 0);   // +2 along y
    check_translator(tx * ty, 1, -0.5, -1.0, 0, "translators add their ideal parts");
    check_translator(ty * tx, 1, -0.5, -1.0, 0, "translators commute");

    // Two planes multiply to a motor; motor_product packs it.
    const Vector<double> p = test_plane(1, 0, 0, 0), q = test_plane(0, 1, 0, -1);
    check_mv_close(ga::to_multivector(ga::motor_product(p, q)), p * q, "motor_product(plane, plane) == plane * plane");
    // Two points too: grades 0 and 2 only.
    const Trivector<double> P = test_point(1, 2, 3), Q = test_point(-1, 0, 2);
    const Multivector<double> PQ = P * Q;
    static_assert(std::is_same<decltype(P * Q), Multivector<double>>::value, "Trivector * Trivector is a Multivector in PGA");
    check(PQ.vector.magnitude().value == 0.0 && PQ.vector.e0 == 0.0 && PQ.trivector.e123 == 0.0
          && PQ.trivector.e012 == 0.0 && PQ.trivector.e013 == 0.0 && PQ.trivector.e023 == 0.0
          && PQ.quadvector.e0123 == 0.0, "point * point has grades 0 and 2 only");
    check_scalar(PQ.scalar, -1.0, "unit points multiply to scalar -1 plus an ideal line");
    check_mv_close(ga::to_multivector(ga::motor_product(P, Q)), PQ, "motor_product(point, point) == point * point");
}

void test_mixed_products()
{
    section("mixed products");

    const Vector<double> v(1, -2, 0.5, 3);
    const Bivector<double> b(2, -1, 0.5, 1, 3, -2);
    const Trivector<double> t(-1, 2, 0.5, 1.5);
    const Quadvector<double> q(2.5);
    const Multivector<double> m = make_mv(1, -2, 3, 0.5, -1, 2, -3, 1, 0.25, -0.5, 4, 1, -1, 2, -2, 0.75);
    const Rotor<double> r(Scalar<double>(0.5), 0.25, -0.75, 1.0);
    const Translator<double> tr(Scalar<double>(1.0), 0.5, -1.0, 2.0);
    const Motor<double> mo(Scalar<double>(0.75), Bivector<double>(1, -2, 0.5, 3, -1, 2), Quadvector<double>(-0.5));

    const Multivector<double> V = ga::detail::promote(v), B = ga::detail::promote(b), Tt = ga::detail::promote(t);
    const Multivector<double> Q = ga::detail::promote(q), R = ga::detail::promote(r), Tr = ga::detail::promote(tr);
    const Multivector<double> Mo = ga::detail::promote(mo);

    // A representative sweep of the 58 pairs against the table, both orders.
    check_mv_close(v * b, V * B, "Vector * Bivector");
    check_mv_close(b * v, B * V, "Bivector * Vector");
    check_mv_close(v * t, V * Tt, "Vector * Trivector");
    check_mv_close(t * v, Tt * V, "Trivector * Vector");
    check_mv_close(v * q, V * Q, "Vector * Quadvector");
    check_mv_close(q * v, Q * V, "Quadvector * Vector");
    check_mv_close(b * b, B * B, "Bivector * Bivector");
    check_mv_close(b * t, B * Tt, "Bivector * Trivector");
    check_mv_close(t * b, Tt * B, "Trivector * Bivector");
    check_mv_close(b * q, B * Q, "Bivector * Quadvector");
    check_mv_close(q * b, Q * B, "Quadvector * Bivector");
    check_mv_close(t * t, Tt * Tt, "Trivector * Trivector");
    check_mv_close(t * q, Tt * Q, "Trivector * Quadvector");
    check_mv_close(q * t, Q * Tt, "Quadvector * Trivector");
    check_mv_close(m * v, m * V, "Multivector * Vector");
    check_mv_close(v * m, V * m, "Vector * Multivector");
    check_mv_close(m * b, m * B, "Multivector * Bivector");
    check_mv_close(t * m, Tt * m, "Trivector * Multivector");
    check_mv_close(m * q, m * Q, "Multivector * Quadvector");
    check_mv_close(q * m, Q * m, "Quadvector * Multivector");
    check_mv_close(r * v, R * V, "Rotor * Vector");
    check_mv_close(v * r, V * R, "Vector * Rotor");
    check_mv_close(tr * t, Tr * Tt, "Translator * Trivector");
    check_mv_close(t * tr, Tt * Tr, "Trivector * Translator");
    check_mv_close(mo * b, Mo * B, "Motor * Bivector");
    check_mv_close(b * mo, B * Mo, "Bivector * Motor");
    check_mv_close(mo * q, Mo * Q, "Motor * Quadvector");
    check_mv_close(m * mo, m * Mo, "Multivector * Motor");
    check_mv_close(r * m, R * m, "Rotor * Multivector");
    check_mv_close(tr * m, Tr * m, "Translator * Multivector");

    // A line squares to a study number: scalar plus e0123, nothing else.
    const Multivector<double> bb = b * b;
    check(bb.vector.e0 == 0.0 && bb.vector.e1 == 0.0 && bb.bivector.e01 == 0.0 && bb.bivector.e12 == 0.0
          && bb.trivector.e123 == 0.0, "a bivector squared has grades 0 and 4 only");
    check_scalar(bb.scalar, -(1 * 1 + 3 * 3 + 2 * 2), "the scalar part of B*B is minus the Euclidean norm squared");
    check_quadvector(bb.quadvector, 2 * (2 * -2 - (-1) * 3 + 0.5 * 1), "the e0123 part of B*B is 2(e01 e23 - e02 e13 + e03 e12)");

    // Higher wedges via the table.
    check_trivector(v ^ b, (V * B).trivector.e012, (V * B).trivector.e013, (V * B).trivector.e023, (V * B).trivector.e123,
                    "Vector ^ Bivector is the grade 3 part");
    check_trivector(b ^ v, (B * V).trivector.e012, (B * V).trivector.e013, (B * V).trivector.e023, (B * V).trivector.e123,
                    "Bivector ^ Vector is the grade 3 part");
    check_quadvector(v ^ t, (V * Tt).quadvector.e0123, "Vector ^ Trivector is the grade 4 part");
    check_quadvector(t ^ v, (Tt * V).quadvector.e0123, "Trivector ^ Vector is the grade 4 part");
    check_quadvector(b ^ b, (B * B).quadvector.e0123, "Bivector ^ Bivector is the grade 4 part");
    check_quadvector(b ^ Bivector<double>(0, 0, 0, 1, 0, 0), (B * ga::detail::promote(Bivector<double>(0, 0, 0, 1, 0, 0))).quadvector.e0123,
                     "Bivector ^ Bivector against a second line");
    check_trivector(ga::wedge_product(v, b), (v ^ b).e012, (v ^ b).e013, (v ^ b).e023, (v ^ b).e123, "wedge_product(v, b) == v ^ b");
    check_trivector(ga::wedge_product(b, v), (b ^ v).e012, (b ^ v).e013, (b ^ v).e023, (b ^ v).e123, "wedge_product(b, v) == b ^ v");
    check_quadvector(ga::wedge_product(v, t), (v ^ t).e0123, "wedge_product(v, t) == v ^ t");
    check_quadvector(ga::wedge_product(t, v), (t ^ v).e0123, "wedge_product(t, v) == t ^ v");
    check_quadvector(ga::wedge_product(b, b), (b ^ b).e0123, "wedge_product(b, b) == b ^ b");

    // The wedge of a plane with a line is the point where they meet.
    const Vector<double> z1 = test_plane(0, 0, 1, -1);
    const Bivector<double> zAxis(0, 0, 0, 1, 0, 0);
    const Trivector<double> hit = z1 ^ zAxis;
    check_close(hit.e012 / hit.e123, -1.0, "z=1 meets the z-axis at z = 1 (e012 = -z)");
    check_close(hit.e013 / hit.e123, 0.0, "... with y = 0");
    check_close(hit.e023 / hit.e123, 0.0, "... with x = 0");
}

void test_mixed_grade_addition()
{
    section("mixed-grade addition");

    const Scalar<double> s(1.0);
    const Vector<double> v(2, 3, 4, 5);
    const Bivector<double> b(6, 7, 8, 9, 10, 11);
    const Trivector<double> t(12, 13, 14, 15);
    const Quadvector<double> q(16.0);
    const Multivector<double> all = make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);

    // A chain through every grade, in an arbitrary order, lands every component.
    check_mv_close(t + s + q + v + b, all, "s + v + b + t + q in any order");

    // Every distinct pair, both orders, against promote-and-add.
    check_mv_close(s + v, make_mv(1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), "Scalar + Vector");
    check_mv_close(v + s, s + v, "Vector + Scalar");
    check_mv_close(s + b, make_mv(1, 0, 0, 0, 0, 6, 7, 8, 9, 10, 11, 0, 0, 0, 0, 0), "Scalar + Bivector");
    check_mv_close(b + s, s + b, "Bivector + Scalar");
    check_mv_close(s + t, make_mv(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 13, 14, 15, 0), "Scalar + Trivector");
    check_mv_close(t + s, s + t, "Trivector + Scalar");
    check_mv_close(s + q, make_mv(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16), "Scalar + Quadvector");
    check_mv_close(q + s, s + q, "Quadvector + Scalar");
    check_mv_close(v + b, make_mv(0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 0, 0, 0, 0), "Vector + Bivector");
    check_mv_close(b + v, v + b, "Bivector + Vector");
    check_mv_close(v + t, make_mv(0, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 12, 13, 14, 15, 0), "Vector + Trivector");
    check_mv_close(t + v, v + t, "Trivector + Vector");
    check_mv_close(v + q, make_mv(0, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16), "Vector + Quadvector");
    check_mv_close(q + v, v + q, "Quadvector + Vector");
    check_mv_close(b + t, make_mv(0, 0, 0, 0, 0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0), "Bivector + Trivector");
    check_mv_close(t + b, b + t, "Trivector + Bivector");
    check_mv_close(b + q, make_mv(0, 0, 0, 0, 0, 6, 7, 8, 9, 10, 11, 0, 0, 0, 0, 16), "Bivector + Quadvector");
    check_mv_close(q + b, b + q, "Quadvector + Bivector");
    check_mv_close(t + q, make_mv(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 13, 14, 15, 16), "Trivector + Quadvector");
    check_mv_close(q + t, t + q, "Quadvector + Trivector");

    // Multivector with each grade, both orders.
    check_mv_close(all + s, make_mv(2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), "Multivector + Scalar");
    check_mv_close(s + all, all + s, "Scalar + Multivector");
    check_mv_close(all + v, make_mv(1, 4, 6, 8, 10, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), "Multivector + Vector");
    check_mv_close(v + all, all + v, "Vector + Multivector");
    check_mv_close(all + b, make_mv(1, 2, 3, 4, 5, 12, 14, 16, 18, 20, 22, 12, 13, 14, 15, 16), "Multivector + Bivector");
    check_mv_close(b + all, all + b, "Bivector + Multivector");
    check_mv_close(all + t, make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 24, 26, 28, 30, 16), "Multivector + Trivector");
    check_mv_close(t + all, all + t, "Trivector + Multivector");
    check_mv_close(all + q, make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 32), "Multivector + Quadvector");
    check_mv_close(q + all, all + q, "Quadvector + Multivector");

    // Every mixed sum is a Multivector, even the ones that would fit a motor.
    static_assert(std::is_same<decltype(s + b), Multivector<double>>::value, "Scalar + Bivector widens to Multivector");
    static_assert(std::is_same<decltype(s + q), Multivector<double>>::value, "Scalar + Quadvector widens to Multivector");
    static_assert(std::is_same<decltype(v + t), Multivector<double>>::value, "Vector + Trivector widens to Multivector");
    static_assert(std::is_same<decltype(all + q), Multivector<double>>::value, "Multivector + Quadvector stays Multivector");

    // The even types widen into Multivector on the right.
    const Rotor<double> r(Scalar<double>(1.0), 2, 3, 4);
    const Translator<double> tr(Scalar<double>(1.0), 2, 3, 4);
    const Motor<double> mo(Scalar<double>(1.0), Bivector<double>(1, 1, 1, 1, 1, 1), Quadvector<double>(1.0));
    check_mv_close(all + r, make_mv(2, 2, 3, 4, 5, 6, 7, 8, 11, 13, 15, 12, 13, 14, 15, 16), "Multivector + Rotor");
    check_mv_close(all + tr, make_mv(2, 2, 3, 4, 5, 8, 10, 12, 9, 10, 11, 12, 13, 14, 15, 16), "Multivector + Translator");
    check_mv_close(all + mo, make_mv(2, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 12, 13, 14, 15, 17), "Multivector + Motor");

    // The named conversions.
    const Motor<double> rm = ga::to_motor(r);
    check_bivector(rm.bivector, 0, 0, 0, 2, 3, 4, "to_motor(Rotor)");
    const Motor<double> tm = ga::to_motor(tr);
    check_bivector(tm.bivector, 2, 3, 4, 0, 0, 0, "to_motor(Translator)");
    const Motor<double> am = ga::to_motor(all);
    check_scalar(am.scalar, 1.0, "to_motor(Multivector) scalar");
    check_bivector(am.bivector, 6, 7, 8, 9, 10, 11, "to_motor(Multivector) bivector");
    check_quadvector(am.quadvector, 16.0, "to_motor(Multivector) quadvector");
    check_rotor(ga::to_rotor(all), 1, 9, 10, 11, "to_rotor(Multivector)");
    check_rotor(ga::to_rotor(mo), 1, 1, 1, 1, "to_rotor(Motor)");
    check_translator(ga::to_translator(all), 1, 6, 7, 8, "to_translator(Multivector)");
    check_translator(ga::to_translator(mo), 1, 1, 1, 1, "to_translator(Motor)");
    check_mv_close(ga::to_multivector(r), make_mv(1, 0, 0, 0, 0, 0, 0, 0, 2, 3, 4, 0, 0, 0, 0, 0), "to_multivector(Rotor)");
    check_mv_close(ga::to_multivector(tr), make_mv(1, 0, 0, 0, 0, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0), "to_multivector(Translator)");
    check_mv_close(ga::to_multivector(mo), make_mv(1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1), "to_multivector(Motor)");
}

void test_mixed_grade_subtraction()
{
    section("mixed-grade subtraction");

    const Scalar<double> s(1.0);
    const Vector<double> v(2, 3, 4, 5);
    const Bivector<double> b(6, 7, 8, 9, 10, 11);
    const Trivector<double> t(12, 13, 14, 15);
    const Quadvector<double> q(16.0);
    const Multivector<double> all = make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);

    // a - b == a + (-b) for every pair, both orders.
    check_mv_close(s - v, s + (-v), "Scalar - Vector");     check_mv_close(v - s, v + (-s), "Vector - Scalar");
    check_mv_close(s - b, s + (-b), "Scalar - Bivector");   check_mv_close(b - s, b + (-s), "Bivector - Scalar");
    check_mv_close(s - t, s + (-t), "Scalar - Trivector");  check_mv_close(t - s, t + (-s), "Trivector - Scalar");
    check_mv_close(s - q, s + (-q), "Scalar - Quadvector"); check_mv_close(q - s, q + (-s), "Quadvector - Scalar");
    check_mv_close(v - b, v + (-b), "Vector - Bivector");   check_mv_close(b - v, b + (-v), "Bivector - Vector");
    check_mv_close(v - t, v + (-t), "Vector - Trivector");  check_mv_close(t - v, t + (-v), "Trivector - Vector");
    check_mv_close(v - q, v + (-q), "Vector - Quadvector"); check_mv_close(q - v, q + (-v), "Quadvector - Vector");
    check_mv_close(b - t, b + (-t), "Bivector - Trivector"); check_mv_close(t - b, t + (-b), "Trivector - Bivector");
    check_mv_close(b - q, b + (-q), "Bivector - Quadvector"); check_mv_close(q - b, q + (-b), "Quadvector - Bivector");
    check_mv_close(t - q, t + (-q), "Trivector - Quadvector"); check_mv_close(q - t, q + (-t), "Quadvector - Trivector");
    check_mv_close(all - s, all + (-s), "Multivector - Scalar");     check_mv_close(s - all, s + (-all), "Scalar - Multivector");
    check_mv_close(all - v, all + (-v), "Multivector - Vector");     check_mv_close(v - all, v + (-all), "Vector - Multivector");
    check_mv_close(all - b, all + (-b), "Multivector - Bivector");   check_mv_close(b - all, b + (-all), "Bivector - Multivector");
    check_mv_close(all - t, all + (-t), "Multivector - Trivector");  check_mv_close(t - all, t + (-all), "Trivector - Multivector");
    check_mv_close(all - q, all + (-q), "Multivector - Quadvector"); check_mv_close(q - all, q + (-all), "Quadvector - Multivector");

    // Subtracting a grade removes it.
    check_mv_close(all - v - b - t - q - s, Multivector<double>(), "removing every grade leaves zero");
    static_assert(std::is_same<decltype(v - q), Multivector<double>>::value, "Vector - Quadvector widens to Multivector");
}

void test_grade_projection()
{
    section("grade projection");

    const Multivector<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
    check_scalar(ga::grade0(m), 1.0, "grade0");
    check_vector(ga::grade1(m), 2, 3, 4, 5, "grade1");
    check_bivector(ga::grade2(m), 6, 7, 8, 9, 10, 11, "grade2");
    check_trivector(ga::grade3(m), 12, 13, 14, 15, "grade3");
    check_quadvector(ga::grade4(m), 16.0, "grade4");
    check_mv_close(ga::grade0(m) + ga::grade1(m) + ga::grade2(m) + ga::grade3(m) + ga::grade4(m), m,
                   "the grades reassemble the multivector");

    const Rotor<double> r(Scalar<double>(1.0), 2, 3, 4);
    check_scalar(ga::grade0(r), 1.0, "grade0(Rotor)");
    check_bivector(ga::grade2(r), 0, 0, 0, 2, 3, 4, "grade2(Rotor) has zero ideal components");

    const Translator<double> t(Scalar<double>(1.0), 2, 3, 4);
    check_scalar(ga::grade0(t), 1.0, "grade0(Translator)");
    check_bivector(ga::grade2(t), 2, 3, 4, 0, 0, 0, "grade2(Translator) has zero Euclidean components");

    const Motor<double> mo(Scalar<double>(1.0), Bivector<double>(2, 3, 4, 5, 6, 7), Quadvector<double>(8.0));
    check_scalar(ga::grade0(mo), 1.0, "grade0(Motor)");
    check_bivector(ga::grade2(mo), 2, 3, 4, 5, 6, 7, "grade2(Motor)");
    check_quadvector(ga::grade4(mo), 8.0, "grade4(Motor)");
}

// ---------------------------------------------------------------------------
// Instantiation sweep
//
// Every check above runs on double and asserts about VALUES. A member of a
// class template is only compiled when it is called, so an entry point nothing
// calls can be syntactically broken and still ship green -- which is exactly
// how Cl3's to_string() once reached main() broken on five of six types. This
// section CALLS every public entry point for float, double and long double so
// a compile failure anywhere in the API stops the build. Add to it whenever
// you add API.
// ---------------------------------------------------------------------------

template <typename T>
void instantiate_every_entry_point()
{
    ga::Scalar<T> s(2);
    ga::Vector<T> v(1, 2, 3, 4);
    ga::Bivector<T> b(1, 2, 3, 4, 5, 6);
    ga::Trivector<T> t(1, 2, 3, 4);
    ga::Quadvector<T> q(5);
    ga::Multivector<T> m(s, v, b, t, q);
    ga::Rotor<T> r(ga::Scalar<T>(1), 0, 0, 0);
    ga::Translator<T> tr(ga::Scalar<T>(1), 0, 0, 0);
    ga::Motor<T> mo(s, b, q);

    // Same-type arithmetic and unary minus.
    (void)(s + s); (void)(s - s); (void)(s * s); (void)(s / s); (void)(-s);
    (void)(v + v); (void)(v - v); (void)(-v);
    (void)(b + b); (void)(b - b); (void)(-b);
    (void)(t + t); (void)(t - t); (void)(-t);
    (void)(q + q); (void)(q - q); (void)(-q); (void)(q * q);
    (void)(m + m); (void)(m - m); (void)(-m);
    (void)(r + r); (void)(r - r); (void)(-r);
    (void)(tr + tr); (void)(tr - tr); (void)(-tr);
    (void)(mo + mo); (void)(mo - mo); (void)(-mo);

    // Scalar mixing, both spellings and both sides.
    (void)(s * v); (void)(s * b); (void)(s * t); (void)(s * q); (void)(s * m); (void)(s * r); (void)(s * tr); (void)(s * mo);
    (void)(v * s); (void)(b * s); (void)(t * s); (void)(q * s); (void)(m * s); (void)(r * s); (void)(tr * s); (void)(mo * s);
    (void)(v / s); (void)(b / s); (void)(t / s); (void)(q / s); (void)(m / s); (void)(r / s); (void)(tr / s); (void)(mo / s);
    (void)(q * T(2)); (void)(q / T(2));
    (void)(T(2) * v); (void)(T(2) * b); (void)(T(2) * t); (void)(T(2) * q); (void)(T(2) * m);
    (void)(T(2) * r); (void)(T(2) * tr); (void)(T(2) * mo);

    // Magnitudes on the blade types.
    (void)v.magnitude(); (void)b.magnitude(); (void)t.magnitude();

    // Compound assignment, every operator on every type.
    { ga::Scalar<T> x = s; x += s; x -= s; x *= s; x /= s; }
    { ga::Vector<T> x = v; x += v; x -= v; x *= T(2); x /= T(2); x *= s; x /= s; }
    { ga::Bivector<T> x = b; x += b; x -= b; x *= T(2); x /= T(2); x *= s; x /= s; }
    { ga::Trivector<T> x = t; x += t; x -= t; x *= T(2); x /= T(2); x *= s; x /= s; }
    { ga::Quadvector<T> x = q; x += q; x -= q; x *= T(2); x /= T(2); x *= s; x /= s; }
    { ga::Multivector<T> x = m; x += m; x -= m; x *= s; x /= s; }
    { ga::Rotor<T> x = r; x += r; x -= r; x *= s; x /= s; }
    { ga::Translator<T> x = tr; x += tr; x -= tr; x *= s; x /= s; }
    { ga::Motor<T> x = mo; x += mo; x -= mo; x *= s; x /= s; }

    // Conversions: implicit widening and explicit narrowing.
    { ga::Multivector<T> x = r; (void)x; }
    { ga::Multivector<T> x = tr; (void)x; }
    { ga::Multivector<T> x = mo; (void)x; }
    { ga::Motor<T> x = r; (void)x; }
    { ga::Motor<T> x = tr; (void)x; }
    { ga::Rotor<T> x(m); (void)x; }
    { ga::Rotor<T> x(mo); (void)x; }
    { ga::Translator<T> x(m); (void)x; }
    { ga::Translator<T> x(mo); (void)x; }
    { ga::Motor<T> x(m); (void)x; }

    // Products: dot, wedge, the general product, and named spellings.
    (void)(v | v); (void)ga::dot_product(v, v);
    (void)(v ^ v); (void)ga::wedge_product(v, v);
    (void)(v ^ b); (void)(b ^ v); (void)(v ^ t); (void)(t ^ v); (void)(b ^ b);
    (void)ga::wedge_product(v, b); (void)ga::wedge_product(b, v); (void)ga::wedge_product(v, t);
    (void)ga::wedge_product(t, v); (void)ga::wedge_product(b, b);
    (void)(v * v); (void)ga::geometric_product(v, v); (void)ga::motor_product(v, v);
    (void)(m * m); (void)ga::geometric_product(m, m);
    (void)(r * r); (void)ga::rotor_product(r, r);
    (void)(tr * tr); (void)ga::translator_product(tr, tr);
    (void)(mo * mo); (void)ga::motor_product(mo, mo);
    (void)ga::motor_product(t, t);

    // Every mixed_products.hpp overload, both orders.
    (void)(v * b); (void)(b * v); (void)(v * t); (void)(t * v); (void)(v * q); (void)(q * v);
    (void)(v * m); (void)(m * v); (void)(v * r); (void)(r * v); (void)(v * tr); (void)(tr * v); (void)(v * mo); (void)(mo * v);
    (void)(b * b); (void)(b * t); (void)(t * b); (void)(b * q); (void)(q * b);
    (void)(b * m); (void)(m * b); (void)(b * r); (void)(r * b); (void)(b * tr); (void)(tr * b); (void)(b * mo); (void)(mo * b);
    (void)(t * t); (void)(t * q); (void)(q * t);
    (void)(t * m); (void)(m * t); (void)(t * r); (void)(r * t); (void)(t * tr); (void)(tr * t); (void)(t * mo); (void)(mo * t);
    (void)(q * m); (void)(m * q); (void)(q * r); (void)(r * q); (void)(q * tr); (void)(tr * q); (void)(q * mo); (void)(mo * q);
    (void)(m * r); (void)(r * m); (void)(m * tr); (void)(tr * m); (void)(m * mo); (void)(mo * m);
    (void)(r * tr); (void)(tr * r); (void)(r * mo); (void)(mo * r); (void)(tr * mo); (void)(mo * tr);

    // Mixed-grade sums and differences, both orders.
    (void)(s + v); (void)(v + s); (void)(s + b); (void)(b + s); (void)(s + t); (void)(t + s); (void)(s + q); (void)(q + s);
    (void)(v + b); (void)(b + v); (void)(v + t); (void)(t + v); (void)(v + q); (void)(q + v);
    (void)(b + t); (void)(t + b); (void)(b + q); (void)(q + b); (void)(t + q); (void)(q + t);
    (void)(m + s); (void)(s + m); (void)(m + v); (void)(v + m); (void)(m + b); (void)(b + m);
    (void)(m + t); (void)(t + m); (void)(m + q); (void)(q + m);
    (void)(s - v); (void)(v - s); (void)(s - b); (void)(b - s); (void)(s - t); (void)(t - s); (void)(s - q); (void)(q - s);
    (void)(v - b); (void)(b - v); (void)(v - t); (void)(t - v); (void)(v - q); (void)(q - v);
    (void)(b - t); (void)(t - b); (void)(b - q); (void)(q - b); (void)(t - q); (void)(q - t);
    (void)(m - s); (void)(s - m); (void)(m - v); (void)(v - m); (void)(m - b); (void)(b - m);
    (void)(m - t); (void)(t - m); (void)(m - q); (void)(q - m);
    (void)(m + r); (void)(m + tr); (void)(m + mo); (void)(mo + r); (void)(mo + tr);

    // Named conversions and grade projection.
    (void)ga::to_motor(r); (void)ga::to_motor(tr); (void)ga::to_motor(m);
    (void)ga::to_rotor(m); (void)ga::to_rotor(mo); (void)ga::to_translator(m); (void)ga::to_translator(mo);
    (void)ga::to_multivector(r); (void)ga::to_multivector(tr); (void)ga::to_multivector(mo);
    (void)ga::grade0(m); (void)ga::grade1(m); (void)ga::grade2(m); (void)ga::grade3(m); (void)ga::grade4(m);
    (void)ga::grade0(r); (void)ga::grade2(r); (void)ga::grade0(tr); (void)ga::grade2(tr);
    (void)ga::grade0(mo); (void)ga::grade2(mo); (void)ga::grade4(mo);

    // Printing.
    (void)s.to_string(); (void)v.to_string(); (void)b.to_string(); (void)t.to_string(); (void)q.to_string();
    (void)m.to_string(); (void)r.to_string(); (void)tr.to_string(); (void)mo.to_string();
}

void test_instantiation_sweep()
{
    section("instantiation sweep");

    instantiate_every_entry_point<float>();
    instantiate_every_entry_point<double>();
    instantiate_every_entry_point<long double>();
    check(true, "every public entry point instantiates for float, double and long double");

    const Scalar<double> s(1.5);
    const Vector<double> v(1, 2, 3, 4);
    const Bivector<double> b(1, 2, 3, 4, 5, 6);
    const Trivector<double> t(1, 2, 3, 4);
    const Quadvector<double> q(5);
    const Multivector<double> m(s, v, b, t, q);
    const Rotor<double> r(s, 1, 2, 3);
    const Translator<double> tr(s, 1, 2, 3);
    const Motor<double> mo(s, b, q);

    check(!s.to_string().empty(), "Scalar::to_string is non-empty");
    check(!v.to_string().empty(), "Vector::to_string is non-empty");
    check(!b.to_string().empty(), "Bivector::to_string is non-empty");
    check(!t.to_string().empty(), "Trivector::to_string is non-empty");
    check(!q.to_string().empty(), "Quadvector::to_string is non-empty");
    check(!m.to_string().empty(), "Multivector::to_string is non-empty");
    check(!r.to_string().empty(), "Rotor::to_string is non-empty");
    check(!tr.to_string().empty(), "Translator::to_string is non-empty");
    check(!mo.to_string().empty(), "Motor::to_string is non-empty");

    // The basis labels must match the component names in docs/pga.md, or the
    // output misleads about which blade is which.
    check(v.to_string().find("e0") != std::string::npos, "Vector labels its e0 component");
    check(b.to_string().find("e01") != std::string::npos, "Bivector labels its e01 component");
    check(b.to_string().find("e23") != std::string::npos, "Bivector labels its e23 component");
    check(t.to_string().find("e012") != std::string::npos, "Trivector labels its e012 component");
    check(q.to_string().find("e0123") != std::string::npos, "Quadvector labels its e0123 component");
    check(r.to_string().find("e12") != std::string::npos, "Rotor labels its e12 component");
    check(tr.to_string().find("e01") != std::string::npos, "Translator labels its e01 component");

    // Round-trip precision, as in the Cl3 suite.
    check(Scalar<double>(1e-17).to_string() != "0.000000",
          "to_string keeps small values instead of flattening them to zero");
    check(Scalar<double>(1e-17).to_string().find("e-17") != std::string::npos,
          "to_string uses scientific notation where it is needed");
}

} // namespace

int main()
{
    test_scalar();
    test_vector();
    test_bivector();
    test_trivector();
    test_quadvector();
    test_multivector();
    test_rotor();
    test_translator();
    test_motor();
    test_compound_assignment();

    test_basis_multiplication_table();
    test_cayley_table_against_permutation_rule();
    test_multivector_product();
    test_dot_and_wedge();
    test_incidence();
    test_even_closure();
    test_mixed_products();
    test_mixed_grade_addition();
    test_mixed_grade_subtraction();
    test_grade_projection();

    test_instantiation_sweep();

    std::cout << g_checks << " checks, " << g_failures << " failed.\n";
    if (g_failures != 0) {
        std::cout << "TESTS FAILED\n";
        return 1;
    }
    std::cout << "All CliffordCore PGA tests passed\n";
    return 0;
}
