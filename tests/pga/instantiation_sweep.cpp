// Cl(3,0,1) tests -- instantiation_sweep: every public entry point called in
// float, double and long double, because a template only fails when it is
// instantiated.
//
// One category of the PGA suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

// ---------------------------------------------------------------------------
// Instantiation sweep
//
// Every other check in this suite runs on double and asserts about VALUES. A
// member of a class template is only compiled when it is called, so an entry
// point nothing calls can be syntactically broken and still ship green --
// which is exactly how Cl3's to_string() once reached main() broken on five of
// six types. This file CALLS every public entry point for float, double and
// long double so a compile failure anywhere in the API stops the build. Add to
// it whenever you add API.
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

    // Norms, on every type.
    (void)ga::norm(s); (void)ga::norm(v); (void)ga::norm(b); (void)ga::norm(t); (void)ga::norm(q);
    (void)ga::norm(m); (void)ga::norm(r); (void)ga::norm(tr); (void)ga::norm(mo);
    (void)ga::squared_norm(s); (void)ga::squared_norm(v); (void)ga::squared_norm(b); (void)ga::squared_norm(t); (void)ga::squared_norm(q);
    (void)ga::squared_norm(m); (void)ga::squared_norm(r); (void)ga::squared_norm(tr); (void)ga::squared_norm(mo);
    (void)ga::ideal_norm(s); (void)ga::ideal_norm(v); (void)ga::ideal_norm(b); (void)ga::ideal_norm(t); (void)ga::ideal_norm(q);
    (void)ga::ideal_norm(m); (void)ga::ideal_norm(r); (void)ga::ideal_norm(tr); (void)ga::ideal_norm(mo);
    (void)ga::normalize(s); (void)ga::normalize(v); (void)ga::normalize(b); (void)ga::normalize(t);
    (void)ga::normalize(m); (void)ga::normalize(r); (void)ga::normalize(tr); (void)ga::normalize(mo);

    // Involutions and inverses.
    (void)ga::reverse(s); (void)ga::reverse(v); (void)ga::reverse(b); (void)ga::reverse(t); (void)ga::reverse(q);
    (void)ga::reverse(m); (void)ga::reverse(r); (void)ga::reverse(tr); (void)ga::reverse(mo);
    (void)ga::involute(s); (void)ga::involute(v); (void)ga::involute(b); (void)ga::involute(t); (void)ga::involute(q);
    (void)ga::involute(m); (void)ga::involute(r); (void)ga::involute(tr); (void)ga::involute(mo);
    (void)ga::conjugate(s); (void)ga::conjugate(v); (void)ga::conjugate(b); (void)ga::conjugate(t); (void)ga::conjugate(q);
    (void)ga::conjugate(m); (void)ga::conjugate(r); (void)ga::conjugate(tr); (void)ga::conjugate(mo);
    (void)ga::inverse(s); (void)ga::inverse(v); (void)ga::inverse(b); (void)ga::inverse(t);
    (void)ga::inverse(m); (void)ga::inverse(r); (void)ga::inverse(tr); (void)ga::inverse(mo);

    // Duality, join and meet.
    (void)ga::dual(s); (void)ga::dual(v); (void)ga::dual(b); (void)ga::dual(t); (void)ga::dual(q); (void)ga::dual(m);
    (void)ga::undual(s); (void)ga::undual(v); (void)ga::undual(b); (void)ga::undual(t); (void)ga::undual(q); (void)ga::undual(m);
    (void)ga::regressive_product(t, t); (void)ga::regressive_product(t, b); (void)ga::regressive_product(b, t); (void)ga::regressive_product(b, b);
    (void)ga::join(t, t); (void)ga::join(t, b); (void)ga::join(b, t); (void)ga::join(b, b);
    (void)ga::meet(v, v); (void)ga::meet(v, b); (void)ga::meet(b, v);

    // Contractions.
    (void)ga::left_contraction(v, v); (void)ga::left_contraction(v, b); (void)ga::left_contraction(v, t); (void)ga::left_contraction(v, q);
    (void)ga::left_contraction(b, b); (void)ga::left_contraction(b, t); (void)ga::left_contraction(b, q);
    (void)ga::left_contraction(t, t); (void)ga::left_contraction(t, q); (void)ga::left_contraction(q, q);
    (void)(v | b); (void)(v | t); (void)(v | q); (void)(b | b); (void)(b | t); (void)(b | q); (void)(t | t); (void)(t | q); (void)(q | q);
    (void)ga::right_contraction(v, v); (void)ga::right_contraction(b, v); (void)ga::right_contraction(t, v); (void)ga::right_contraction(q, v);
    (void)ga::right_contraction(b, b); (void)ga::right_contraction(t, b); (void)ga::right_contraction(q, b);
    (void)ga::right_contraction(t, t); (void)ga::right_contraction(q, t); (void)ga::right_contraction(q, q);
    (void)ga::scalar_product(v, v); (void)ga::scalar_product(b, b); (void)ga::scalar_product(m, m); (void)ga::scalar_product(v, b);
    (void)ga::scalar_product(q, q); (void)ga::scalar_product(s, s); (void)ga::grade0(s);

    // Comparison, on every type.
    (void)(s == s); (void)(v == v); (void)(b == b); (void)(t == t); (void)(q == q); (void)(m == m); (void)(r == r); (void)(tr == tr); (void)(mo == mo);
    (void)(s != s); (void)(v != v); (void)(b != b); (void)(t != t); (void)(q != q); (void)(m != m); (void)(r != r); (void)(tr != tr); (void)(mo != mo);
    (void)ga::approx_equal(s, s); (void)ga::approx_equal(v, v); (void)ga::approx_equal(b, b); (void)ga::approx_equal(t, t);
    (void)ga::approx_equal(q, q); (void)ga::approx_equal(m, m); (void)ga::approx_equal(r, r); (void)ga::approx_equal(tr, tr); (void)ga::approx_equal(mo, mo);
    (void)ga::approx_equal(v, v, T(1e-3));

    // Primitives and the accessors.
    (void)ga::point(T(1), T(2), T(3)); (void)ga::ideal_point(T(1), T(2), T(3)); (void)ga::origin<T>();
    (void)ga::plane(T(1), T(2), T(3), T(4)); (void)ga::ideal_plane<T>(); (void)ga::line_through_origin(T(1), T(2), T(3));
    (void)ga::line_through_points(t, t); (void)ga::line_from_planes(v, v); (void)ga::plane_through_points(t, t, t);
    (void)ga::position(t); (void)ga::direction(t); (void)ga::direction(b); (void)ga::moment(b);
    (void)ga::normal(v); (void)ga::offset(v);
    (void)ga::is_ideal(v); (void)ga::is_ideal(b); (void)ga::is_ideal(t);

    // exp, log and the sandwich.
    (void)ga::exp(b); (void)ga::log(mo); (void)ga::log(r); (void)ga::log(tr);
    (void)ga::sandwich(v, mo); (void)ga::sandwich(b, mo); (void)ga::sandwich(t, mo); (void)ga::sandwich(m, mo);
    (void)ga::sandwich(v, r); (void)ga::sandwich(b, r); (void)ga::sandwich(t, r); (void)ga::sandwich(m, r);
    (void)ga::sandwich(v, tr); (void)ga::sandwich(b, tr); (void)ga::sandwich(t, tr); (void)ga::sandwich(m, tr);
    (void)ga::rotate(v, r); (void)ga::rotate(b, r); (void)ga::rotate(t, r); (void)ga::rotate(m, r);
    (void)ga::translate(v, tr); (void)ga::translate(b, tr); (void)ga::translate(t, tr); (void)ga::translate(m, tr);
    (void)ga::transform(v, mo); (void)ga::transform(b, mo); (void)ga::transform(t, mo); (void)ga::transform(m, mo);

    // Versor construction.
    (void)ga::identity_rotor<T>(); (void)ga::identity_translator<T>(); (void)ga::identity_motor<T>();
    (void)ga::translator(T(1), T(2), T(3)); (void)ga::translator_between(t, t);
    (void)ga::rotor_from_axis_angle(T(0), T(0), T(1), T(1));
    (void)ga::motor_from_line_angle(b, T(1)); (void)ga::screw(b, T(1), T(2));
    (void)ga::motor_between(v, v); (void)ga::motor_from_rotor_translator(r, tr);
    (void)ga::slerp(mo, mo, T(0.5)); (void)ga::slerp(r, r, T(0.5)); (void)ga::slerp(tr, tr, T(0.5));

    // Geometry.
    (void)ga::reflect(v, v); (void)ga::reflect(b, v); (void)ga::reflect(t, v);
    (void)ga::project(t, v); (void)ga::project(t, b); (void)ga::project(b, v);
    (void)ga::distance(t, t); (void)ga::distance(t, v); (void)ga::distance(t, b);
    (void)ga::angle(v, v); (void)ga::angle(b, b);

    // Printing.
    (void)s.to_string(); (void)v.to_string(); (void)b.to_string(); (void)t.to_string(); (void)q.to_string();
    (void)m.to_string(); (void)r.to_string(); (void)tr.to_string(); (void)mo.to_string();
    {
        std::ostringstream os;
        os << s << v << b << t << q << m << r << tr << mo;
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

void run_instantiation_sweep()
{
    test_instantiation_sweep();
}

} // namespace cctest
