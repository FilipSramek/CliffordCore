// Cl(3,0) tests -- instantiation_sweep: every public entry point called in
// float, double and long double, because a template only fails when it is
// instantiated.
//
// One category of the Cl3 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

// ---------------------------------------------------------------------------
// Instantiation sweep
//
// Every other check in this suite runs on double and asserts about VALUES.
// That leaves a blind spot: a member of a class template is only compiled when
// it is called, so an entry point nothing calls can be syntactically broken and
// still ship green. That is exactly how to_string() reached main() broken on
// five of the six types while this suite reported "779 checks, 0 failed".
//
// This file closes it. It asserts little -- its job is to CALL every public
// entry point, for float, double and long double, so a compile failure
// anywhere in the API stops the build.
// ---------------------------------------------------------------------------

template <typename T>
void instantiate_every_entry_point()
{
    ga::Scalar<T> s(2);
    ga::Vector<T> v(1, 2, 3);
    ga::Bivector<T> b(1, 2, 3);
    ga::Trivector<T> t(4);
    ga::Multivector<T> m(s, v, b, t);
    ga::Rotor<T> r(ga::Scalar<T>(1), ga::Bivector<T>(0, 0, 0));

    // Same-type arithmetic and unary minus.
    (void)(s + s); (void)(s - s); (void)(s * s); (void)(s / s); (void)(-s);
    (void)(v + v); (void)(v - v); (void)(-v);
    (void)(b + b); (void)(b - b); (void)(-b);
    (void)(t + t); (void)(t - t); (void)(-t); (void)(t * t); (void)(t / t);
    (void)(m + m); (void)(m - m); (void)(-m);
    (void)(r + r); (void)(r - r);

    // Scalar mixing, both spellings and both sides.
    (void)(s * v); (void)(s * b); (void)(s * t); (void)(s * m); (void)(s * r);
    (void)(v * s); (void)(b * s); (void)(t * s); (void)(m * s); (void)(r * s);
    (void)(v / s); (void)(b / s); (void)(t / s); (void)(m / s); (void)(r / s);
    (void)(t * T(2)); (void)(t / T(2));
    (void)(T(2) * v); (void)(T(2) * b); (void)(T(2) * t); (void)(T(2) * m); (void)(T(2) * r);

    // Compound assignment, every operator on every type.
    { ga::Scalar<T> q = s; q += s; q -= s; q *= s; q /= s; }
    { ga::Vector<T> q = v; q += v; q -= v; q *= T(2); q /= T(2); q *= s; q /= s; }
    { ga::Bivector<T> q = b; q += b; q -= b; q *= T(2); q /= T(2); q *= s; q /= s; }
    { ga::Trivector<T> q = t; q += t; q -= t; q *= T(2); q /= T(2); q *= s; q /= s; }
    { ga::Multivector<T> q = m; q += m; q -= m; q *= s; q /= s; }
    { ga::Rotor<T> q = r; q += r; q -= r; q *= s; q /= s; }

    // Products.
    (void)(v | v); (void)ga::dot_product(v, v);
    (void)(v ^ v); (void)ga::wedge_product(v, v);
    (void)(v ^ b); (void)(b ^ v); (void)ga::wedge_product(v, b);
    (void)(v * v); (void)(m * m); (void)(r * r);
    (void)ga::geometric_product(v, v);
    (void)ga::geometric_product(m, m);
    (void)ga::rotor_product(v, v);
    (void)ga::rotor_product(r, r);

    // Every mixed_products.hpp overload, in both orders.
    (void)(v * b); (void)(b * v);
    (void)(b * b);
    (void)(v * t); (void)(t * v);
    (void)(b * t); (void)(t * b);
    (void)(r * v); (void)(v * r);
    (void)(m * v); (void)(v * m);
    (void)(m * b); (void)(b * m);
    (void)(m * t); (void)(t * m);
    (void)(m * r); (void)(r * m);

    // Mixed-grade sums and differences, both orders.
    (void)(s + v); (void)(v + s); (void)(s + b); (void)(b + s); (void)(s + t); (void)(t + s);
    (void)(v + b); (void)(b + v); (void)(v + t); (void)(t + v); (void)(b + t); (void)(t + b);
    (void)(m + s); (void)(s + m); (void)(m + v); (void)(v + m);
    (void)(m + b); (void)(b + m); (void)(m + t); (void)(t + m);
    (void)(s - v); (void)(v - s); (void)(s - b); (void)(b - s); (void)(s - t); (void)(t - s);
    (void)(v - b); (void)(b - v); (void)(v - t); (void)(t - v); (void)(b - t); (void)(t - b);
    (void)(m - s); (void)(s - m); (void)(m - v); (void)(v - m);
    (void)(m - b); (void)(b - m); (void)(m - t); (void)(t - m);

    // Magnitudes, involutions, structure.
    (void)ga::norm(s); (void)ga::norm(v); (void)ga::norm(b);
    (void)ga::norm(t); (void)ga::norm(m); (void)ga::norm(r);
    (void)ga::squared_norm(v); (void)ga::squared_norm(b);
    (void)ga::squared_norm(t); (void)ga::squared_norm(m);
    (void)ga::squared_norm(r);
    (void)ga::normalize(v); (void)ga::normalize(b);
    (void)ga::normalize(t); (void)ga::normalize(m);
    (void)ga::normalize(r);
    (void)ga::reverse(v); (void)ga::reverse(b);
    (void)ga::reverse(t); (void)ga::reverse(m);
    (void)ga::reverse(r);
    (void)ga::inverse(s); (void)ga::inverse(v);
    (void)ga::inverse(b); (void)ga::inverse(t);
    (void)ga::inverse(m); (void)ga::inverse(r);
    (void)ga::conjugate(m);
    (void)ga::dual(s); (void)ga::dual(v);
    (void)ga::dual(b); (void)ga::dual(t);
    (void)ga::grade0(m); (void)ga::grade1(m);
    (void)ga::grade2(m); (void)ga::grade3(m);
    (void)ga::grade0(r); (void)ga::grade2(r);
    (void)ga::to_rotor(m); (void)ga::to_multivector(r);
    (void)ga::rotor_sum(s, b);
    (void)v.magnitude(); (void)b.magnitude(); (void)t.magnitude();

    // Rotations.
    (void)ga::exp(b);
    (void)ga::log(r);
    (void)ga::sandwich(v, r); (void)ga::rotate(v, r);
    (void)ga::identity_rotor<T>();
    (void)ga::rotor_from_axis_angle(v, T(1));
    (void)ga::rotor_between(v, v);
    (void)ga::slerp(r, r, T(0.5));

    // Printing -- the entry points that shipped broken because nothing called them.
    (void)s.to_string(); (void)v.to_string(); (void)b.to_string();
    (void)t.to_string(); (void)m.to_string(); (void)r.to_string();

    // The four operator* pairs added for the bivector sandwich.
    (void)(r * b); (void)(b * r); (void)(r * t); (void)(t * r);

    // Involutions, on every type.
    (void)ga::involute(s); (void)ga::involute(v);
    (void)ga::involute(b); (void)ga::involute(t);
    (void)ga::involute(m); (void)ga::involute(r);
    (void)ga::conjugate(s); (void)ga::conjugate(v);
    (void)ga::conjugate(b); (void)ga::conjugate(t);
    (void)ga::conjugate(m); (void)ga::conjugate(r);
    (void)ga::reverse(s);
    (void)ga::squared_norm(s); (void)ga::normalize(s);

    // Contractions and the scalar product.
    (void)ga::left_contraction(v, v); (void)ga::left_contraction(v, b);
    (void)ga::left_contraction(v, t); (void)ga::left_contraction(b, b);
    (void)ga::left_contraction(b, t); (void)ga::left_contraction(t, t);
    (void)ga::right_contraction(v, v); (void)ga::right_contraction(b, v);
    (void)ga::right_contraction(t, v); (void)ga::right_contraction(b, b);
    (void)ga::right_contraction(t, b); (void)ga::right_contraction(t, t);
    (void)(v | b); (void)(v | t); (void)(b | b); (void)(b | t);
    (void)ga::scalar_product(v, v); (void)ga::scalar_product(b, b);
    (void)ga::scalar_product(t, t); (void)ga::scalar_product(s, s); (void)ga::grade0(s);

    // Geometry.
    (void)ga::reflect(v, v); (void)ga::reflect(v, b);
    (void)ga::project(v, v); (void)ga::project(v, b);
    (void)ga::reject(v, v); (void)ga::reject(v, b);

    // Sandwich, on every grade it accepts.
    (void)ga::sandwich(b, r); (void)ga::sandwich(t, r);
    (void)ga::sandwich(m, r);
    (void)ga::rotate(b, r); (void)ga::rotate(t, r);
    (void)ga::rotate(m, r);

    // Comparison, exact and tolerant.
    (void)(s == s); (void)(v == v); (void)(b == b);
    (void)(t == t); (void)(m == m); (void)(r == r);
    (void)(s != s); (void)(v != v); (void)(b != b);
    (void)(t != t); (void)(m != m); (void)(r != r);
    (void)ga::approx_equal(s, s); (void)ga::approx_equal(v, v);
    (void)ga::approx_equal(b, b); (void)ga::approx_equal(t, t);
    (void)ga::approx_equal(m, m); (void)ga::approx_equal(r, r);

    // Stream insertion, for all six.
    {
        std::ostringstream os;
        os << s << v << b << t << m << r;
        (void)os;
    }
}

void test_instantiation_sweep()
{
    section("instantiation sweep");

    instantiate_every_entry_point<float>();
    instantiate_every_entry_point<double>();
    instantiate_every_entry_point<long double>();
    check(true, "every public entry point instantiates for float, double and long double");

    const Scalar<double> s(1.5);
    const Vector<double> v(1, 2, 3);
    const Bivector<double> b(1, 2, 3);
    const Trivector<double> t(4);
    const Multivector<double> m(s, v, b, t);
    const Rotor<double> r(s, b);

    check(!s.to_string().empty(), "Scalar::to_string is non-empty");
    check(!v.to_string().empty(), "Vector::to_string is non-empty");
    check(!b.to_string().empty(), "Bivector::to_string is non-empty");
    check(!t.to_string().empty(), "Trivector::to_string is non-empty");
    check(!m.to_string().empty(), "Multivector::to_string is non-empty");
    check(!r.to_string().empty(), "Rotor::to_string is non-empty");

    // The basis labels must match the component names in docs/conventions.md,
    // or the output misleads about which grade is which.
    check(v.to_string().find("e1") != std::string::npos, "Vector labels its e1 component");
    check(b.to_string().find("e12") != std::string::npos, "Bivector labels its e12 component");
    check(t.to_string().find("e123") != std::string::npos, "Trivector labels its e123 component");

    // Round-trip precision. std::to_string renders this as "0.000000", which is
    // useless for the residuals this library actually produces.
    check(Scalar<double>(1e-17).to_string() != "0.000000",
          "to_string keeps small values instead of flattening them to zero");
    check(Scalar<double>(1e-17).to_string().find("e-17") != std::string::npos,
          "to_string uses scientific notation where it is needed");
}

} // namespace

void run_instantiation_sweep()
{
    test_instantiation_sweep();
}

} // namespace cctest
