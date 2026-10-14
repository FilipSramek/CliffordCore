// Cl(3,0) tests -- arithmetic: mixed-grade sums and differences, compound
// assignment, and raw scalars on the left.
//
// One category of the Cl3 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

void test_mixed_grade_addition()
{
    section("mixed-grade addition");

    const Scalar<double> s(1.0);
    const Vector<double> v1(1, 0, 0), v2(0, 1, 0), v3(0, 0, 1), v4(1, 1, 1);
    const Bivector<double> b1(1, 0, 0), b2(0, 1, 0), b3(0, 0, 1);
    const Trivector<double> t(5.0);

    // A chain across every grade, in one expression.
    const Multivector<double> m = s + v1 + v2 + v3 + v4 + b1 + b2 + b3 + t;
    check_mv(m, 1, 2, 2, 2, 1, 1, 1, 5, "nine-term mixed chain");

    // Addition is commutative, so the order of the chain must not matter.
    const Multivector<double> reordered = t + b3 + b2 + b1 + v4 + v3 + v2 + v1 + s;
    check_close(mv_difference(m, reordered), 0.0, "chain order does not matter");

    // Every mixed pair lands in the right slot.
    check_mv(s + v1, 1, 1, 0, 0, 0, 0, 0, 0, "Scalar + Vector");
    check_mv(v1 + s, 1, 1, 0, 0, 0, 0, 0, 0, "Vector + Scalar");
    check_mv(s + b1, 1, 0, 0, 0, 1, 0, 0, 0, "Scalar + Bivector");
    check_mv(s + t, 1, 0, 0, 0, 0, 0, 0, 5, "Scalar + Trivector");
    check_mv(v1 + b1, 0, 1, 0, 0, 1, 0, 0, 0, "Vector + Bivector");
    check_mv(v1 + t, 0, 1, 0, 0, 0, 0, 0, 5, "Vector + Trivector");
    check_mv(b1 + t, 0, 0, 0, 0, 1, 0, 0, 5, "Bivector + Trivector");
    check_mv(m + v1, 1, 3, 2, 2, 1, 1, 1, 5, "Multivector + Vector");
    check_mv(v1 + m, 1, 3, 2, 2, 1, 1, 1, 5, "Vector + Multivector");

    // Same-grade sums must not widen.
    static_assert(std::is_same<decltype(v1 + v2), Vector<double>>::value,
                  "Vector + Vector stays a Vector");
    static_assert(std::is_same<decltype(s + b1), Multivector<double>>::value,
                  "Scalar + Bivector widens to Multivector");

    // rotor_sum packs the same numbers as s + b, and to_rotor narrows.
    const Rotor<double> viaSum = ga::rotor_sum(s, b1);
    check_scalar(viaSum.scalar, 1.0, "rotor_sum scalar");
    check_bivector(viaSum.bivector, 1.0, 0.0, 0.0, "rotor_sum bivector");
    const Rotor<double> viaNarrow = ga::to_rotor(s + b1);
    check_scalar(viaNarrow.scalar, 1.0, "to_rotor scalar");
    check_bivector(viaNarrow.bivector, 1.0, 0.0, 0.0, "to_rotor bivector");
    check_mv(ga::to_multivector(viaSum), 1, 0, 0, 0, 1, 0, 0, 0, "to_multivector(rotor)");
}

void test_mixed_grade_subtraction()
{
    section("mixed-grade subtraction");

    const Scalar<double> s(2.0);
    const Vector<double> v(1, 2, 3);
    const Bivector<double> b(4, 5, 6);
    const Trivector<double> t(7.0);
    const Multivector<double> m = make_mv(1, 1, 2, 3, 1, 2, 3, 4);

    check_mv(s - v, 2, -1, -2, -3, 0, 0, 0, 0, "Scalar - Vector");
    check_mv(v - s, -2, 1, 2, 3, 0, 0, 0, 0, "Vector - Scalar");
    check_mv(s - b, 2, 0, 0, 0, -4, -5, -6, 0, "Scalar - Bivector");
    check_mv(s - t, 2, 0, 0, 0, 0, 0, 0, -7, "Scalar - Trivector");
    check_mv(v - b, 0, 1, 2, 3, -4, -5, -6, 0, "Vector - Bivector");
    check_mv(b - v, 0, -1, -2, -3, 4, 5, 6, 0, "Bivector - Vector");
    check_mv(v - t, 0, 1, 2, 3, 0, 0, 0, -7, "Vector - Trivector");
    check_mv(b - t, 0, 0, 0, 0, 4, 5, 6, -7, "Bivector - Trivector");
    check_mv(m - v, 1, 0, 0, 0, 1, 2, 3, 4, "Multivector - Vector");
    check_mv(v - m, -1, 0, 0, 0, -1, -2, -3, -4, "Vector - Multivector");

    // a - b must equal a + (-b), and a - a must vanish.
    check_close(mv_difference(v - b, v + (-b)), 0.0, "a - b == a + (-b)");
    check_close(mv_difference(m - m, Multivector<double>()), 0.0, "m - m == 0");

    // Same-grade differences must not widen.
    static_assert(std::is_same<decltype(v - Vector<double>()), Vector<double>>::value,
                  "Vector - Vector stays a Vector");
    static_assert(std::is_same<decltype(s - v), Multivector<double>>::value,
                  "Scalar - Vector widens to Multivector");
}

void test_compound_assignment()
{
    section("compound assignment");

    Vector<double> v(1, 0, 0);
    v += Vector<double>(0, 1, 0);
    check_vector(v, 1.0, 1.0, 0.0, "Vector +=");
    v -= Vector<double>(1, 0, 0);
    check_vector(v, 0.0, 1.0, 0.0, "Vector -=");
    v *= 3.0;
    check_vector(v, 0.0, 3.0, 0.0, "Vector *= raw");
    v /= Scalar<double>(3.0);
    check_vector(v, 0.0, 1.0, 0.0, "Vector /= Scalar");

    Bivector<double> b(1, 2, 3);
    b += Bivector<double>(1, 1, 1);
    check_bivector(b, 2.0, 3.0, 4.0, "Bivector +=");
    b *= Scalar<double>(2.0);
    check_bivector(b, 4.0, 6.0, 8.0, "Bivector *= Scalar");

    Trivector<double> t(4.0);
    t += Trivector<double>(2.0);
    check_trivector(t, 6.0, "Trivector +=");
    t /= 2.0;
    check_trivector(t, 3.0, "Trivector /= raw");

    Scalar<double> s(10.0);
    s -= Scalar<double>(4.0);
    check_scalar(s, 6.0, "Scalar -=");
    s *= Scalar<double>(2.0);
    check_scalar(s, 12.0, "Scalar *=");

    Multivector<double> m = make_mv(1, 1, 1, 1, 1, 1, 1, 1);
    m += make_mv(1, 1, 1, 1, 1, 1, 1, 1);
    check_mv(m, 2, 2, 2, 2, 2, 2, 2, 2, "Multivector +=");
    m *= Scalar<double>(0.5);
    check_mv(m, 1, 1, 1, 1, 1, 1, 1, 1, "Multivector *= Scalar");

    Rotor<double> r(Scalar<double>(1.0), Bivector<double>(2.0, 0.0, 0.0));
    r += Rotor<double>(Scalar<double>(1.0), Bivector<double>(2.0, 0.0, 0.0));
    check_scalar(r.scalar, 2.0, "Rotor += scalar part");
    check_bivector(r.bivector, 4.0, 0.0, 0.0, "Rotor += bivector part");
}

void test_scalar_on_the_left()
{
    section("raw scalar on the left");

    // These used to work only with the number on the right.
    check_vector(3.0 * Vector<double>(1, 2, 3), 3.0, 6.0, 9.0, "raw * Vector");
    check_bivector(3.0 * Bivector<double>(1, 2, 3), 3.0, 6.0, 9.0, "raw * Bivector");
    check_trivector(3.0 * Trivector<double>(2.0), 6.0, "raw * Trivector");

    // Both spellings must agree.
    const Vector<double> v(1, 2, 3);
    const Vector<double> left = 2.0 * v;
    const Vector<double> right = v * 2.0;
    check_vector(left, right.x, right.y, right.z, "2*v == v*2");

    check_mv(2.0 * make_mv(1, 1, 1, 1, 1, 1, 1, 1), 2, 2, 2, 2, 2, 2, 2, 2, "raw * Multivector");
    const Rotor<double> r = 2.0 * Rotor<double>(Scalar<double>(1.0), Bivector<double>(1, 0, 0));
    check_scalar(r.scalar, 2.0, "raw * Rotor");
}

} // namespace

void run_arithmetic()
{
    test_mixed_grade_addition();
    test_mixed_grade_subtraction();
    test_compound_assignment();
    test_scalar_on_the_left();
}

} // namespace cctest
