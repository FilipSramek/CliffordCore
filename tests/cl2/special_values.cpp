// Cl(2,0) tests -- special_values: NaN, zero and the edges of the input domain.
//
// The same three questions as tests/cl3/special_values.cpp, which explains
// them (docs/safety-audit.md, items A1, A3 and A4): NaN must stay contagious,
// zero inputs are pinned at today's behaviour, and results must stay accurate
// at both edges of the documented input domain.
//
// One category of the Cl2 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

// True if any component is NaN or infinite.
template<typename T> bool any_nonfinite(T x) { return !std::isfinite(x); }
template<typename T> bool any_nonfinite(const Scalar<T>& s) { return !std::isfinite(s.value); }
template<typename T> bool any_nonfinite(const Vector<T>& v) { return !std::isfinite(v.x) || !std::isfinite(v.y); }
template<typename T> bool any_nonfinite(const Bivector<T>& b) { return !std::isfinite(b.xy); }
template<typename T> bool any_nonfinite(const Rotor<T>& r)
{
    return any_nonfinite(r.scalar) || any_nonfinite(r.bivector);
}
template<typename T> bool any_nonfinite(const Multivector<T>& m)
{
    return any_nonfinite(m.scalar) || any_nonfinite(m.vector) || any_nonfinite(m.bivector);
}

// The input domain, as documented in docs/safety-audit.md.
template<typename T> T domain_max();
template<> float domain_max<float>() { return std::ldexp(1.0f, 30); }
template<> double domain_max<double>() { return std::ldexp(1.0, 250); }
template<> long double domain_max<long double>() { return std::ldexp(1.0L, 250); }

// ---------------------------------------------------------------------------
// NaN contagion
// ---------------------------------------------------------------------------

template<typename T>
void test_nan_contagion()
{
    static std::string label;
    label = std::string("NaN contagion, ") + type_name<T>();
    section(label.c_str());

    const T nan = std::numeric_limits<T>::quiet_NaN();
    const Vector<T> w(T(0.3), T(-0.5));
    const Bivector<T> b2(T(0.4));
    const Rotor<T> r2 = ga::rotor_from_angle(T(0.9));
    const Multivector<T> m2(Scalar<T>(T(1)), w, b2);

    for (int i = 0; i < 2; ++i) {
        Vector<T> v(T(1), T(2));
        (i == 0 ? v.x : v.y) = nan;
        check(any_nonfinite(ga::norm(v)), "norm(Vector)");
        check(any_nonfinite(ga::squared_norm(v)), "squared_norm(Vector)");
        check(any_nonfinite(ga::normalize(v)), "normalize(Vector)");
        check(any_nonfinite(ga::inverse(v)), "inverse(Vector)");
        check(any_nonfinite(ga::reverse(v)), "reverse(Vector)");
        check(any_nonfinite(ga::involute(v)), "involute(Vector)");
        check(any_nonfinite(ga::conjugate(v)), "conjugate(Vector)");
        check(any_nonfinite(ga::dual(v)), "dual(Vector)");
        check(any_nonfinite(ga::perp(v)), "perp(Vector)");
        check(any_nonfinite(v ^ w), "Vector ^ Vector");
        check(any_nonfinite(v | w), "Vector | Vector");
        check(any_nonfinite(v * w), "Vector * Vector");
        check(any_nonfinite(v * b2), "Vector * Bivector");
        check(any_nonfinite(ga::rotate(v, r2)), "rotate(Vector, Rotor)");
        check(any_nonfinite(ga::reflect(v, w)), "reflect(v, n)");
        check(any_nonfinite(ga::reflect(w, v)), "reflect(v, n) with NaN in n");
        check(any_nonfinite(ga::project(v, w)), "project(v, u)");
        check(any_nonfinite(ga::project(w, v)), "project(v, u) with NaN in u");
        check(any_nonfinite(ga::reject(v, w)), "reject(v, u)");
        check(any_nonfinite(ga::signed_area(v, w)), "signed_area");
        check(any_nonfinite(ga::angle(v)), "angle(Vector)");
        check(any_nonfinite(ga::angle_between(v, w)), "angle_between(NaN, w)");
        check(any_nonfinite(ga::angle_between(w, v)), "angle_between(w, NaN)");
        check(any_nonfinite(ga::rotor_between(v, w)), "rotor_between(NaN, w)");
        check(any_nonfinite(ga::rotor_between(w, v)), "rotor_between(w, NaN)");
    }

    const Bivector<T> b(nan);
    check(any_nonfinite(ga::norm(b)), "norm(Bivector)");
    check(any_nonfinite(ga::normalize(b)), "normalize(Bivector)");
    check(any_nonfinite(ga::inverse(b)), "inverse(Bivector)");
    check(any_nonfinite(ga::dual(b)), "dual(Bivector)");
    check(any_nonfinite(ga::exp(b)), "exp(Bivector)");
    check(any_nonfinite(b * b2), "Bivector * Bivector");

    for (int i = 0; i < 2; ++i) {
        Rotor<T> r = r2;
        (i == 0 ? r.scalar.value : r.bivector.xy) = nan;
        check(any_nonfinite(ga::norm(r)), "norm(Rotor)");
        check(any_nonfinite(ga::normalize(r)), "normalize(Rotor)");
        check(any_nonfinite(ga::inverse(r)), "inverse(Rotor)");
        check(any_nonfinite(ga::reverse(r)), "reverse(Rotor)");
        check(any_nonfinite(ga::log(r)), "log(Rotor)");
        check(any_nonfinite(ga::rotor_angle(r)), "rotor_angle(Rotor)");
        check(any_nonfinite(r * r2), "Rotor * Rotor");
        check(any_nonfinite(ga::rotate(w, r)), "rotate(v, NaN rotor)");
        check(any_nonfinite(ga::rotate(m2, r)), "rotate(Multivector, NaN rotor)");
        check(any_nonfinite(ga::slerp(r, r2, T(0.5))), "slerp(NaN, r, 0.5)");
        check(any_nonfinite(ga::slerp(r2, r, T(0.5))), "slerp(r, NaN, 0.5)");
    }

    for (int i = 0; i < 4; ++i) {
        Multivector<T> m = m2;
        T* const slot[4] = {&m.scalar.value, &m.vector.x, &m.vector.y, &m.bivector.xy};
        *slot[i] = nan;
        check(any_nonfinite(ga::norm(m)), "norm(Multivector)");
        check(any_nonfinite(ga::normalize(m)), "normalize(Multivector)");
        check(any_nonfinite(ga::inverse(m)), "inverse(Multivector)");
        check(any_nonfinite(ga::dual(m)), "dual(Multivector)");
        check(any_nonfinite(m * m2), "Multivector * Multivector");
        check(any_nonfinite(ga::rotate(m, r2)), "rotate(Multivector, Rotor)");
    }

    check(any_nonfinite(ga::rotor_from_angle(nan)), "rotor_from_angle(NaN)");
    check(any_nonfinite(ga::slerp(ga::identity_rotor<T>(), r2, nan)), "slerp(a, b, NaN)");
    check(any_nonfinite(ga::exp(Bivector<T>(std::numeric_limits<T>::infinity()))), "exp(infinite Bivector)");
}

// ---------------------------------------------------------------------------
// Zero inputs: today's behaviour, pinned
// ---------------------------------------------------------------------------

// The zeros below reach the library through a volatile, so the optimizer
// cannot fold them into a constant division by zero (MSVC /O2 warns C4723).
// See the same helper in tests/cl3/special_values.cpp.
double opaque_zero()
{
    volatile double zero = 0.0;
    return zero;
}

void test_zero_inputs()
{
    section("zero inputs (today's behaviour; see docs/safety-audit.md A1, A3)");

    const double z = opaque_zero();
    const Vector<double> zero(z, z);
    const Vector<double> e1(1.0, 0.0);
    // Copy-initialised: with two Type(name) arguments, direct initialisation
    // would parse as a function declaration.
    const Rotor<double> zeroRotor = Rotor<double>(Scalar<double>(z), Bivector<double>(z));

    // A1: a degenerate input returns a value that looks like a real answer.
    check(ga::rotor_between(zero, e1) == ga::identity_rotor<double>(), "rotor_between(zero, v) is the identity rotor");
    check(ga::rotor_between(e1, zero) == ga::identity_rotor<double>(), "rotor_between(v, zero) is the identity rotor");
    check(ga::angle(zero) == 0.0, "angle(zero vector) is 0");
    check(ga::angle_between(zero, e1) == 0.0, "angle_between(zero, v) is 0");
    check(ga::rotor_angle(zeroRotor) == 0.0, "rotor_angle(zero rotor) is 0");
    check(ga::log(zeroRotor) == Bivector<double>(0.0), "log(zero rotor) is the zero bivector, log of the identity");

    // A3: single-grade inverses divide by zero, the Multivector one does not.
    check(std::isinf(ga::inverse(Scalar<double>(z)).value), "inverse(Scalar 0) is infinite");
    check(std::isnan(ga::inverse(zero).x), "inverse(Vector 0) is NaN");
    check(std::isnan(ga::inverse(Bivector<double>(z)).xy), "inverse(Bivector 0) is NaN");
    check(any_nonfinite(ga::inverse(zeroRotor)), "inverse(Rotor 0) is NaN");
    check(ga::inverse(Multivector<double>()) == Multivector<double>(), "inverse(Multivector 0) is zero");
    check(any_nonfinite(ga::project(e1, zero)), "project onto a zero vector is NaN");
    check(any_nonfinite(ga::reflect(e1, zero)), "reflect in a zero normal is NaN");

    check(ga::exp(Bivector<double>(0.0)) == ga::identity_rotor<double>(), "exp(0) is the identity rotor");
}

// ---------------------------------------------------------------------------
// The edges of the input domain
// ---------------------------------------------------------------------------

template<typename T>
void test_domain_edge(T s)
{
    check_close_t(ga::norm(Vector<T>(T(3) * s, T(4) * s)).value, T(5) * s, "norm(Vector) at the edge");
    const Vector<T> unit = ga::normalize(Vector<T>(T(3) * s, T(4) * s));
    check_close_t(unit.x, T(3) / T(5), "normalize(Vector) at the edge");
    check_close_t(ga::inverse(Vector<T>(s, T(0))).x * s, T(1), "inverse(Vector) at the edge");
    check_close_t(ga::inverse(Bivector<T>(s)).xy * s, T(-1), "inverse(Bivector) at the edge");
    check_close_t(ga::angle(Vector<T>(s, s)), std::atan2(T(1), T(1)), "angle at the edge");

    const Rotor<T> between = ga::rotor_between(Vector<T>(s, T(0)), Vector<T>(T(0), s));
    const Vector<T> turned = ga::rotate(Vector<T>(s, T(0)), between);
    check_close_t(turned.y / s, T(1), "rotor_between at the edge: e1 goes to e2");
    check_close_t(turned.x / s + T(1), T(1), "rotor_between at the edge: no e1 left");
}

template<typename T>
void test_domain_edges()
{
    static std::string label;
    label = std::string("input domain edges, ") + type_name<T>();
    section(label.c_str());

    test_domain_edge<T>(domain_max<T>());
    test_domain_edge<T>(T(1) / domain_max<T>());
}

void test_special_values_in_every_type()
{
    test_nan_contagion<float>();
    test_nan_contagion<double>();
    test_nan_contagion<long double>();
    test_domain_edges<float>();
    test_domain_edges<double>();
    test_domain_edges<long double>();
}

} // namespace

void run_special_values()
{
    test_special_values_in_every_type();
    test_zero_inputs();
}

} // namespace cctest
