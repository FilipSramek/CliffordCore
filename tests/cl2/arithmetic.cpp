// Cl(2,0) tests -- arithmetic: scaling by a Scalar or a raw value on either
// side, compound assignment, and mixed-grade sums and differences.
//
// One category of the Cl2 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

template<typename A, typename B>
using sum_t = decltype(std::declval<A>() + std::declval<B>());

template<typename A, typename B>
using difference_t = decltype(std::declval<A>() - std::declval<B>());

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

} // namespace

void run_arithmetic()
{
    test_scalar_products();
    test_compound_assignment();
    test_mixed_grade_addition();
    test_mixed_grade_subtraction();
}

} // namespace cctest
