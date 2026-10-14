// Cl(2,0) tests -- instantiation_sweep: every public entry point called in
// float, double and long double, because a template only fails when it is
// instantiated.
//
// One category of the Cl2 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

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

void run_instantiation_sweep()
{
    test_instantiation_sweep();
}

} // namespace cctest
