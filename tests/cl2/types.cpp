// Cl(2,0) tests -- types: the five types on their own: construction, same-type
// arithmetic, the bivector as the pseudoscalar, and the conversions between
// Rotor and Multivector.
//
// One category of the Cl2 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

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

} // namespace

void run_types()
{
    test_construction();
    test_same_type_arithmetic();
    test_bivector_is_the_pseudoscalar();
    test_conversions();
}

} // namespace cctest
