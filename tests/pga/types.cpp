// Cl(3,0,1) tests -- types: the nine types on their own: construction and
// same-type arithmetic, plus the compile-time checks on constexpr evaluation,
// widening and narrowing conversions, and the geometric aliases.
//
// One category of the PGA suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

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

} // namespace

void run_types()
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
}

} // namespace cctest
