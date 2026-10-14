// Cl(3,0) tests -- types: the six types on their own: construction and
// same-type arithmetic, plus the compile-time (constexpr) checks.
//
// One category of the Cl3 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

// ---------------------------------------------------------------------------
// Compile-time checks
//
// Everything in the library is constexpr, so these must evaluate during
// compilation. If constexpr-ness regresses, the build fails rather than the run.
// std::sqrt is not portably constexpr, so norm/magnitude stay out of here.
// ---------------------------------------------------------------------------

constexpr Vector<double> kCtA(1.0, 2.0, 3.0);
constexpr Vector<double> kCtB(4.0, 5.0, 6.0);
static_assert((kCtA | kCtB).value == 32.0, "dot product must be constexpr");
static_assert((kCtA ^ kCtB).xy == -3.0, "wedge product must be constexpr");
static_assert(Scalar<double>(2.0).value == 2.0, "Scalar construction must be constexpr");
static_assert(Vector<double>().x == 0.0, "default Vector must be constexpr and zeroed");

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

    // Cross-type multiplication: Scalar on the left, a graded type on the right.
    const Scalar<double> two(2.0);
    check_vector(two * Vector<double>(1.0, 2.0, 3.0), 2.0, 4.0, 6.0, "Scalar * Vector");
    check_bivector(two * Bivector<double>(1.0, 2.0, 3.0), 2.0, 4.0, 6.0, "Scalar * Bivector");
    check_trivector(two * Trivector<double>(5.0), 10.0, "Scalar * Trivector");

    // Composite types: the scalar must distribute over every grade.
    const Multivector<double> mv(
        Scalar<double>(1.0), Vector<double>(2.0, 3.0, 4.0),
        Bivector<double>(5.0, 6.0, 7.0), Trivector<double>(8.0));
    const Multivector<double> scaledMv = two * mv;
    check_scalar(scaledMv.scalar, 2.0, "Scalar * Multivector scalar");
    check_vector(scaledMv.vector, 4.0, 6.0, 8.0, "Scalar * Multivector vector");
    check_bivector(scaledMv.bivector, 10.0, 12.0, 14.0, "Scalar * Multivector bivector");
    check_trivector(scaledMv.trivector, 16.0, "Scalar * Multivector trivector");

    const Rotor<double> rot(Scalar<double>(1.0), Bivector<double>(2.0, 3.0, 4.0));
    const Rotor<double> scaledRot = two * rot;
    check_scalar(scaledRot.scalar, 2.0, "Scalar * Rotor scalar");
    check_bivector(scaledRot.bivector, 4.0, 6.0, 8.0, "Scalar * Rotor bivector");
}

void test_vector3()
{
    section("Vector");

    check_vector(Vector<double>(), 0.0, 0.0, 0.0, "default ctor zeroes");
    check_vector(Vector<double>(1.0, 2.0, 3.0), 1.0, 2.0, 3.0, "component ctor");

    const Vector<double> a(1.0, 2.0, 3.0);
    const Vector<double> b(4.0, 5.0, 6.0);

    check_vector(a + b, 5.0, 7.0, 9.0, "operator+");
    check_vector(a - b, -3.0, -3.0, -3.0, "operator-");
    check_vector(-a, -1.0, -2.0, -3.0, "unary operator-");
    check_vector(a * Scalar<double>(2.0), 2.0, 4.0, 6.0, "operator* (Scalar)");
    check_vector(a / Scalar<double>(2.0), 0.5, 1.0, 1.5, "operator/ (Scalar)");

    // Scalar multiplication agrees across the two spellings.
    const Vector<double> lhs = Scalar<double>(3.0) * a;
    const Vector<double> rhs = a * Scalar<double>(3.0);
    check_vector(lhs, rhs.x, rhs.y, rhs.z, "Scalar*Vector == Vector*Scalar");
}

void test_bivector3()
{
    section("Bivector");

    check_bivector(Bivector<double>(), 0.0, 0.0, 0.0, "default ctor zeroes");
    check_bivector(Bivector<double>(1.0, 2.0, 3.0), 1.0, 2.0, 3.0, "component ctor");

    const Bivector<double> a(1.0, 2.0, 3.0);
    const Bivector<double> b(4.0, 5.0, 6.0);

    check_bivector(a + b, 5.0, 7.0, 9.0, "operator+");
    check_bivector(a - b, -3.0, -3.0, -3.0, "operator-");
    check_bivector(-a, -1.0, -2.0, -3.0, "unary operator-");
    check_bivector(a * Scalar<double>(2.0), 2.0, 4.0, 6.0, "operator* (Scalar)");
    check_bivector(a / Scalar<double>(2.0), 0.5, 1.0, 1.5, "operator/ (Scalar)");

    check_scalar(a.magnitude(), std::sqrt(14.0), "magnitude");
    check_scalar(Bivector<double>().magnitude(), 0.0, "magnitude of zero bivector");

    // magnitude() must agree with the free norm() function.
    check_close(a.magnitude().value, ga::norm(a).value, "magnitude == norm");
}

void test_trivector3()
{
    section("Trivector");

    check_trivector(Trivector<double>(), 0.0, "default ctor zeroes");
    check_trivector(Trivector<double>(7.0), 7.0, "value ctor");

    const Trivector<double> a(6.0);
    const Trivector<double> b(3.0);

    check_trivector(a + b, 9.0, "operator+");
    check_trivector(a - b, 3.0, "operator-");
    check_trivector(-a, -6.0, "unary operator-");

    // The geometric product of two pseudoscalars is a SCALAR, not a trivector:
    // e123 squares to -1, so (6 e123)(3 e123) = -18. Division is a * inverse(b),
    // where the two minus signs cancel and leave a plain 6/3.
    static_assert(std::is_same<decltype(a * b), Scalar<double>>::value,
                  "Trivector * Trivector produces a Scalar");
    static_assert(std::is_same<decltype(a / b), Scalar<double>>::value,
                  "Trivector / Trivector produces a Scalar");
    check_scalar(a * b, -18.0, "operator* (Trivector) is -ab");
    check_scalar(a / b, 2.0, "operator/ (Trivector)");

    // The defining case, and the one the old component-wise version got wrong.
    check_scalar(Trivector<double>(1.0) * Trivector<double>(1.0), -1.0, "e123 * e123 == -1");

    // Division must undo multiplication: (a/b) scaled back by b returns a.
    check_trivector((a / b) * b, a.e123, "(a/b) * b == a");

    check_trivector(a * Scalar<double>(2.0), 12.0, "operator* (Scalar)");
    check_trivector(a / Scalar<double>(2.0), 3.0, "operator/ (Scalar)");

    // Raw numeric overloads exist so `t * 2.0` stays unambiguous now that both
    // Trivector and Scalar are constructible from one T. They must agree with
    // the Scalar spelling.
    check_trivector(a * 2.0, 12.0, "operator* (raw T)");
    check_trivector(a / 2.0, 3.0, "operator/ (raw T)");
    check_close((a * 2.0).e123, (a * Scalar<double>(2.0)).e123, "raw T == Scalar overload");

    check_scalar(a.magnitude(), 6.0, "magnitude");
    check_scalar(Trivector<double>(-6.0).magnitude(), 6.0, "magnitude is non-negative");
}

void test_multivector3()
{
    section("Multivector");

    const Multivector<double> zero;
    check_scalar(zero.scalar, 0.0, "default ctor scalar");
    check_vector(zero.vector, 0.0, 0.0, 0.0, "default ctor vector");
    check_bivector(zero.bivector, 0.0, 0.0, 0.0, "default ctor bivector");
    check_trivector(zero.trivector, 0.0, "default ctor trivector");

    const Multivector<double> a(
        Scalar<double>(1.0), Vector<double>(2.0, 3.0, 4.0),
        Bivector<double>(5.0, 6.0, 7.0), Trivector<double>(8.0));
    const Multivector<double> b(
        Scalar<double>(10.0), Vector<double>(20.0, 30.0, 40.0),
        Bivector<double>(50.0, 60.0, 70.0), Trivector<double>(80.0));

    const Multivector<double> sum = a + b;
    check_scalar(sum.scalar, 11.0, "operator+ scalar");
    check_vector(sum.vector, 22.0, 33.0, 44.0, "operator+ vector");
    check_bivector(sum.bivector, 55.0, 66.0, 77.0, "operator+ bivector");
    check_trivector(sum.trivector, 88.0, "operator+ trivector");

    const Multivector<double> diff = b - a;
    check_scalar(diff.scalar, 9.0, "operator- scalar");
    check_vector(diff.vector, 18.0, 27.0, 36.0, "operator- vector");
    check_bivector(diff.bivector, 45.0, 54.0, 63.0, "operator- bivector");
    check_trivector(diff.trivector, 72.0, "operator- trivector");

    const Multivector<double> negated = -a;
    check_scalar(negated.scalar, -1.0, "unary operator- scalar");
    check_vector(negated.vector, -2.0, -3.0, -4.0, "unary operator- vector");
    check_bivector(negated.bivector, -5.0, -6.0, -7.0, "unary operator- bivector");
    check_trivector(negated.trivector, -8.0, "unary operator- trivector");

    const Multivector<double> scaled = a * Scalar<double>(2.0);
    check_scalar(scaled.scalar, 2.0, "operator* (Scalar) scalar");
    check_vector(scaled.vector, 4.0, 6.0, 8.0, "operator* (Scalar) vector");
    check_bivector(scaled.bivector, 10.0, 12.0, 14.0, "operator* (Scalar) bivector");
    check_trivector(scaled.trivector, 16.0, "operator* (Scalar) trivector");

    const Multivector<double> halved = a / Scalar<double>(2.0);
    check_scalar(halved.scalar, 0.5, "operator/ (Scalar) scalar");
    check_vector(halved.vector, 1.0, 1.5, 2.0, "operator/ (Scalar) vector");
    check_bivector(halved.bivector, 2.5, 3.0, 3.5, "operator/ (Scalar) bivector");
    check_trivector(halved.trivector, 4.0, "operator/ (Scalar) trivector");

    // Scaling agrees whichever side the scalar sits on.
    const Multivector<double> fromLeft = Scalar<double>(2.0) * a;
    check_close(fromLeft.scalar.value, scaled.scalar.value, "Scalar*MV == MV*Scalar scalar");
    check_vector(fromLeft.vector, scaled.vector.x, scaled.vector.y, scaled.vector.z,
                 "Scalar*MV == MV*Scalar vector");
    check_trivector(fromLeft.trivector, scaled.trivector.e123, "Scalar*MV == MV*Scalar trivector");

    // Negation is subtraction from zero.
    const Multivector<double> viaSubtraction = Multivector<double>() - a;
    check_scalar(viaSubtraction.scalar, negated.scalar.value, "-a == 0 - a scalar");
    check_trivector(viaSubtraction.trivector, negated.trivector.e123, "-a == 0 - a trivector");
}

void test_rotor3()
{
    section("Rotor");

    const Rotor<double> zero;
    check_scalar(zero.scalar, 0.0, "default ctor scalar");
    check_bivector(zero.bivector, 0.0, 0.0, 0.0, "default ctor bivector");

    const Rotor<double> a(Scalar<double>(1.0), Bivector<double>(2.0, 3.0, 4.0));
    const Rotor<double> b(Scalar<double>(10.0), Bivector<double>(20.0, 30.0, 40.0));

    const Rotor<double> sum = a + b;
    check_scalar(sum.scalar, 11.0, "operator+ scalar");
    check_bivector(sum.bivector, 22.0, 33.0, 44.0, "operator+ bivector");

    const Rotor<double> diff = b - a;
    check_scalar(diff.scalar, 9.0, "operator- scalar");
    check_bivector(diff.bivector, 18.0, 27.0, 36.0, "operator- bivector");

    const Rotor<double> scaled = a * Scalar<double>(2.0);
    check_scalar(scaled.scalar, 2.0, "operator* (Scalar) scalar");
    check_bivector(scaled.bivector, 4.0, 6.0, 8.0, "operator* (Scalar) bivector");
}

} // namespace

void run_types()
{
    test_scalar();
    test_vector3();
    test_bivector3();
    test_trivector3();
    test_multivector3();
    test_rotor3();
}

} // namespace cctest
