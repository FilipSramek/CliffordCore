// Cl(3,0) tests -- special_values: NaN, zero and the edges of the input domain.
//
// Every other category tests ordinary values. This one tests what a
// safety-critical caller has to know about the inputs that are not ordinary
// (docs/safety-audit.md, items A1, A3 and A4):
//
//   - NaN is contagious: a NaN in any component an operation reads gives a
//     non-finite result, never a plausible finite one. That property is safe
//     to keep under any future error-handling policy, so it is asserted.
//   - Zero inputs: what each degenerate case returns TODAY. Several of these
//     are findings, not designs -- an identity rotor for a zero axis is
//     indistinguishable from a real answer. They are pinned so that the
//     policy change the audit recommends has to update them deliberately.
//   - The documented input domain: every component within [2^-30, 2^30] in
//     float and [2^-250, 2^250] in double and long double. Results at both
//     edges must still be accurate.
//
// One category of the Cl3 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

// True if any component is NaN or infinite.
template<typename T> bool any_nonfinite(const Scalar<T>& s) { return !std::isfinite(s.value); }
template<typename T> bool any_nonfinite(const Vector<T>& v)
{
    return !std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z);
}
template<typename T> bool any_nonfinite(const Bivector<T>& b)
{
    return !std::isfinite(b.xy) || !std::isfinite(b.xz) || !std::isfinite(b.yz);
}
template<typename T> bool any_nonfinite(const Trivector<T>& t) { return !std::isfinite(t.e123); }
template<typename T> bool any_nonfinite(const Rotor<T>& r)
{
    return any_nonfinite(r.scalar) || any_nonfinite(r.bivector);
}
template<typename T> bool any_nonfinite(const Multivector<T>& m)
{
    return any_nonfinite(m.scalar) || any_nonfinite(m.vector)
        || any_nonfinite(m.bivector) || any_nonfinite(m.trivector);
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
    const Vector<T> w(T(0.3), T(-0.5), T(0.7));
    const Bivector<T> b2(T(0.2), T(0.4), T(-0.6));
    const Rotor<T> r2 = ga::rotor_from_axis_angle(Vector<T>(T(1), T(2), T(3)), T(0.9));
    const Multivector<T> m2(Scalar<T>(T(1)), w, b2, Trivector<T>(T(0.5)));

    // A NaN in each component of each operand in turn.
    for (int i = 0; i < 3; ++i) {
        Vector<T> v(T(1), T(2), T(3));
        (i == 0 ? v.x : i == 1 ? v.y : v.z) = nan;
        check(any_nonfinite(ga::norm(v)), "norm(Vector)");
        check(any_nonfinite(ga::squared_norm(v)), "squared_norm(Vector)");
        check(any_nonfinite(ga::normalize(v)), "normalize(Vector)");
        check(any_nonfinite(ga::inverse(v)), "inverse(Vector)");
        check(any_nonfinite(ga::reverse(v)), "reverse(Vector)");
        check(any_nonfinite(ga::involute(v)), "involute(Vector)");
        check(any_nonfinite(ga::conjugate(v)), "conjugate(Vector)");
        check(any_nonfinite(ga::dual(v)), "dual(Vector)");
        check(any_nonfinite(v ^ w), "Vector ^ Vector");
        check(any_nonfinite(v | w), "Vector | Vector");
        check(any_nonfinite(v * w), "Vector * Vector");
        check(any_nonfinite(v ^ b2), "Vector ^ Bivector");
        check(any_nonfinite(ga::left_contraction(v, b2)), "left_contraction(Vector, Bivector)");
        check(any_nonfinite(ga::right_contraction(b2, v)), "right_contraction(Bivector, Vector)");
        check(any_nonfinite(ga::rotate(v, r2)), "rotate(Vector, Rotor)");
        check(any_nonfinite(ga::reflect(v, w)), "reflect(v, n)");
        check(any_nonfinite(ga::reflect(w, v)), "reflect(v, n) with NaN in n");
        check(any_nonfinite(ga::project(v, w)), "project(v, u)");
        check(any_nonfinite(ga::project(w, v)), "project(v, u) with NaN in u");
        check(any_nonfinite(ga::reject(v, w)), "reject(v, u)");
        check(any_nonfinite(ga::project(v, b2)), "project(v, plane)");
        check(any_nonfinite(ga::reject(v, b2)), "reject(v, plane)");
        check(any_nonfinite(ga::rotor_between(v, w)), "rotor_between(NaN, w)");
        check(any_nonfinite(ga::rotor_between(w, v)), "rotor_between(w, NaN)");
        check(any_nonfinite(ga::rotor_from_axis_angle(v, T(0.5))), "rotor_from_axis_angle(NaN axis)");
    }

    for (int i = 0; i < 3; ++i) {
        Bivector<T> b(T(0.1), T(0.2), T(0.3));
        (i == 0 ? b.xy : i == 1 ? b.xz : b.yz) = nan;
        check(any_nonfinite(ga::norm(b)), "norm(Bivector)");
        check(any_nonfinite(ga::normalize(b)), "normalize(Bivector)");
        check(any_nonfinite(ga::inverse(b)), "inverse(Bivector)");
        check(any_nonfinite(ga::dual(b)), "dual(Bivector)");
        check(any_nonfinite(ga::exp(b)), "exp(Bivector)");
        check(any_nonfinite(b * b2), "Bivector * Bivector");
        check(any_nonfinite(ga::rotate(b, r2)), "rotate(Bivector, Rotor)");
        check(any_nonfinite(ga::project(w, b)), "project(v, NaN plane)");
        check(any_nonfinite(ga::reflect(w, b)), "reflect(v, NaN plane)");
    }

    for (int i = 0; i < 4; ++i) {
        Rotor<T> r = r2;
        (i == 0 ? r.scalar.value : i == 1 ? r.bivector.xy : i == 2 ? r.bivector.xz : r.bivector.yz) = nan;
        check(any_nonfinite(ga::norm(r)), "norm(Rotor)");
        check(any_nonfinite(ga::normalize(r)), "normalize(Rotor)");
        check(any_nonfinite(ga::inverse(r)), "inverse(Rotor)");
        check(any_nonfinite(ga::reverse(r)), "reverse(Rotor)");
        check(any_nonfinite(ga::log(r)), "log(Rotor)");
        check(any_nonfinite(r * r2), "Rotor * Rotor");
        check(any_nonfinite(ga::rotate(w, r)), "rotate(v, NaN rotor)");
        check(any_nonfinite(ga::rotate(m2, r)), "rotate(Multivector, NaN rotor)");
        check(any_nonfinite(ga::slerp(r, r2, T(0.5))), "slerp(NaN, r, 0.5)");
        check(any_nonfinite(ga::slerp(r2, r, T(0.5))), "slerp(r, NaN, 0.5)");
    }

    for (int i = 0; i < 8; ++i) {
        Multivector<T> m = m2;
        T* const slot[8] = {&m.scalar.value, &m.vector.x, &m.vector.y, &m.vector.z,
                            &m.bivector.xy, &m.bivector.xz, &m.bivector.yz, &m.trivector.e123};
        *slot[i] = nan;
        check(any_nonfinite(ga::norm(m)), "norm(Multivector)");
        check(any_nonfinite(ga::normalize(m)), "normalize(Multivector)");
        check(any_nonfinite(ga::inverse(m)), "inverse(Multivector)");
        check(any_nonfinite(m * m2), "Multivector * Multivector");
        check(any_nonfinite(ga::rotate(m, r2)), "rotate(Multivector, Rotor)");
    }

    const Trivector<T> t(nan);
    check(any_nonfinite(ga::inverse(t)), "inverse(Trivector)");
    check(any_nonfinite(ga::normalize(t)), "normalize(Trivector)");
    check(any_nonfinite(ga::dual(t)), "dual(Trivector)");
    check(any_nonfinite(t * t), "Trivector * Trivector");

    // A NaN parameter rather than a NaN operand.
    check(any_nonfinite(ga::rotor_from_axis_angle(w, nan)), "rotor_from_axis_angle(axis, NaN angle)");
    check(any_nonfinite(ga::slerp(ga::identity_rotor<T>(), r2, nan)), "slerp(a, b, NaN)");
    check(any_nonfinite(ga::exp(Bivector<T>(std::numeric_limits<T>::infinity(), T(0), T(0)))),
          "exp(infinite Bivector)");
}

// ---------------------------------------------------------------------------
// Zero inputs: today's behaviour, pinned
// ---------------------------------------------------------------------------

// The zeros below reach the library through a volatile, so the optimizer
// cannot fold them into a constant division by zero. MSVC at /O2 otherwise
// warns C4723 at the division inside the library -- at exactly the spot
// these tests exist to pin.
double opaque_zero()
{
    volatile double zero = 0.0;
    return zero;
}

void test_zero_inputs()
{
    section("zero inputs (today's behaviour; see docs/safety-audit.md A1, A3)");

    const double z = opaque_zero();
    const Vector<double> zero(z, z, z);
    const Vector<double> e1(1.0, 0.0, 0.0);
    const Bivector<double> zeroBivector(z, z, z);
    const Rotor<double> zeroRotor(Scalar<double>(z), zeroBivector);

    // A1: a degenerate input returns a value that looks like a real answer.
    check(ga::rotor_from_axis_angle(zero, 1.0) == ga::identity_rotor<double>(),
          "rotor_from_axis_angle(zero axis) is the identity rotor");
    check(ga::rotor_between(zero, e1) == ga::identity_rotor<double>(),
          "rotor_between(zero, v) is the identity rotor");
    check(ga::rotor_between(e1, zero) == ga::identity_rotor<double>(),
          "rotor_between(v, zero) is the identity rotor");
    check(ga::log(zeroRotor) == Bivector<double>(), "log(zero rotor) is the zero bivector, log of the identity");

    // A3: single-grade inverses divide by zero, the Multivector one does not.
    check(std::isinf(ga::inverse(Scalar<double>(z)).value), "inverse(Scalar 0) is infinite");
    check(std::isnan(ga::inverse(zero).x), "inverse(Vector 0) is NaN");
    check(std::isnan(ga::inverse(zeroBivector).xy), "inverse(Bivector 0) is NaN");
    check(std::isnan(ga::inverse(Trivector<double>(z)).e123), "inverse(Trivector 0) is NaN");
    check(any_nonfinite(ga::inverse(zeroRotor)), "inverse(Rotor 0) is NaN");
    check(ga::inverse(Multivector<double>()) == Multivector<double>(), "inverse(Multivector 0) is zero");
    // ...and so do the operations built on inverse.
    check(any_nonfinite(ga::project(e1, zero)), "project onto a zero vector is NaN");
    check(any_nonfinite(ga::reflect(e1, zero)), "reflect in a zero normal is NaN");

    // Correct, not degenerate: exp of the zero bivector really is the identity.
    check(ga::exp(Bivector<double>()) == ga::identity_rotor<double>(), "exp(0) is the identity rotor");
}

// ---------------------------------------------------------------------------
// The edges of the input domain
// ---------------------------------------------------------------------------

template<typename T>
void test_domain_edge(T s)
{
    const T third = T(1) / T(3);

    check_close_t(ga::norm(Vector<T>(s, T(2) * s, T(2) * s)).value, T(3) * s, "norm(Vector) at the edge");
    const Vector<T> unit = ga::normalize(Vector<T>(s, T(2) * s, T(2) * s));
    check_close_t(unit.x, third, "normalize(Vector) at the edge");
    check_close_t(unit.z, T(2) * third, "normalize(Vector) at the edge, z");
    check_close_t(ga::inverse(Vector<T>(s, T(0), T(0))).x * s, T(1), "inverse(Vector) at the edge");
    check_close_t(ga::inverse(Bivector<T>(s, T(0), T(0))).xy * s, T(-1), "inverse(Bivector) at the edge");

    // A quarter turn about z, built from an axis of size s. The angle comes
    // from atan2 in T: the double constant kPi is not accurate to long double.
    const T halfPi = std::atan2(T(1), T(0));
    const Rotor<T> quarter = ga::rotor_from_axis_angle(Vector<T>(T(0), T(0), s), halfPi);
    const Vector<T> turned = ga::rotate(Vector<T>(s, T(0), T(0)), quarter);
    check_close_t(turned.y / s, T(1), "rotate at the edge: e1 goes to e2");
    check_close_t(turned.x / s + T(1), T(1), "rotate at the edge: no e1 left");

    const Rotor<T> between = ga::rotor_between(Vector<T>(s, T(0), T(0)), Vector<T>(T(0), s, T(0)));
    check_close_t(ga::rotate(Vector<T>(T(1), T(0), T(0)), between).y, T(1), "rotor_between at the edge");
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
