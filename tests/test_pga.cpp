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

    test_instantiation_sweep();

    std::cout << g_checks << " checks, " << g_failures << " failed.\n";
    if (g_failures != 0) {
        std::cout << "TESTS FAILED\n";
        return 1;
    }
    std::cout << "All CliffordCore PGA tests passed\n";
    return 0;
}
