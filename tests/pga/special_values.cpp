// Cl(3,0,1) tests -- special_values: NaN, zero and ideal inputs, and the edges
// of the input domain.
//
// The same three questions as tests/cl3/special_values.cpp, which explains
// them (docs/safety-audit.md, items A1, A3 and A4), with one twist the
// degenerate metric adds: norm(), angle() and the inner product ignore every
// e0 component by design, so a NaN in an ideal component is invisible to them.
// The sweep below therefore poisons only Euclidean components for those, and
// a separate test pins that blind spot.
//
// One category of the PGA suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

// True if any component is NaN or infinite.
template<typename T> bool any_nonfinite(T x) { return !std::isfinite(x); }
template<typename T> bool any_nonfinite(const Scalar<T>& s) { return !std::isfinite(s.value); }
template<typename T> bool any_nonfinite(const Vector<T>& v)
{
    return !std::isfinite(v.e0) || !std::isfinite(v.e1) || !std::isfinite(v.e2) || !std::isfinite(v.e3);
}
template<typename T> bool any_nonfinite(const Bivector<T>& b)
{
    return !std::isfinite(b.e01) || !std::isfinite(b.e02) || !std::isfinite(b.e03)
        || !std::isfinite(b.e12) || !std::isfinite(b.e13) || !std::isfinite(b.e23);
}
template<typename T> bool any_nonfinite(const Trivector<T>& t)
{
    return !std::isfinite(t.e012) || !std::isfinite(t.e013) || !std::isfinite(t.e023) || !std::isfinite(t.e123);
}
template<typename T> bool any_nonfinite(const Quadvector<T>& q) { return !std::isfinite(q.e0123); }
template<typename T> bool any_nonfinite(const Multivector<T>& m)
{
    return any_nonfinite(m.scalar) || any_nonfinite(m.vector) || any_nonfinite(m.bivector)
        || any_nonfinite(m.trivector) || any_nonfinite(m.quadvector);
}
template<typename T> bool any_nonfinite(const Rotor<T>& r)
{
    return any_nonfinite(r.scalar) || !std::isfinite(r.e12) || !std::isfinite(r.e13) || !std::isfinite(r.e23);
}
template<typename T> bool any_nonfinite(const Translator<T>& t)
{
    return any_nonfinite(t.scalar) || !std::isfinite(t.e01) || !std::isfinite(t.e02) || !std::isfinite(t.e03);
}
template<typename T> bool any_nonfinite(const Motor<T>& m)
{
    return any_nonfinite(m.scalar) || any_nonfinite(m.bivector) || any_nonfinite(m.quadvector);
}
template<typename T> bool any_nonfinite(const std::array<T, 3>& a)
{
    return !std::isfinite(a[0]) || !std::isfinite(a[1]) || !std::isfinite(a[2]);
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
    const Vector<T> w = ga::plane(T(0.3), T(-0.5), T(0.7), T(1.5));
    const Trivector<T> q = ga::point(T(-1), T(0.5), T(2));
    const Bivector<T> l2 = ga::line_through_points(ga::point(T(1), T(1), T(0)), q);
    const Motor<T> m2 = ga::motor_from_rotor_translator(ga::rotor_from_axis_angle(T(1), T(2), T(3), T(0.9)),
                                                        ga::translator(T(0.4), T(-0.2), T(0.1)));

    // Planes: e0 is ideal, e1..e3 Euclidean.
    for (int i = 0; i < 4; ++i) {
        Vector<T> v = ga::plane(T(1), T(2), T(2), T(-3));
        T* const slot[4] = {&v.e0, &v.e1, &v.e2, &v.e3};
        *slot[i] = nan;
        const bool euclidean = i > 0;
        if (euclidean) {
            // e0 squares to zero, so none of these three reads it.
            check(any_nonfinite(ga::norm(v)), "norm(plane)");
            check(any_nonfinite(ga::angle(v, w)), "angle(plane, plane)");
            check(any_nonfinite(v | w), "plane | plane");
        }
        check(any_nonfinite(ga::normalize(v)), "normalize(plane)");
        check(any_nonfinite(ga::inverse(v)), "inverse(plane)");
        check(any_nonfinite(ga::reverse(v)), "reverse(plane)");
        check(any_nonfinite(ga::dual(v)), "dual(plane)");
        check(any_nonfinite(v ^ w), "plane ^ plane");
        check(any_nonfinite(v * w), "plane * plane");
        check(any_nonfinite(ga::meet(v, l2)), "meet(plane, line)");
        check(any_nonfinite(ga::sandwich(v, m2)), "sandwich(plane, motor)");
        check(any_nonfinite(ga::reflect(v, w)), "reflect(plane, plane)");
        check(any_nonfinite(ga::reflect(q, v)), "reflect(point, NaN plane)");
        check(any_nonfinite(ga::project(q, v)), "project(point, NaN plane)");
        check(any_nonfinite(ga::project(l2, v)), "project(line, NaN plane)");
        check(any_nonfinite(ga::distance(q, v)), "distance(point, NaN plane)");
        check(any_nonfinite(ga::motor_between(v, w)), "motor_between(NaN plane, plane)");
    }

    // Lines: e01..e03 ideal, e12..e23 Euclidean.
    for (int i = 0; i < 6; ++i) {
        Bivector<T> b = l2;
        T* const slot[6] = {&b.e01, &b.e02, &b.e03, &b.e12, &b.e13, &b.e23};
        *slot[i] = nan;
        const bool euclidean = i >= 3;
        if (euclidean) {
            check(any_nonfinite(ga::norm(b)), "norm(line)");
            check(any_nonfinite(ga::angle(b, l2)), "angle(line, line)");
        }
        check(any_nonfinite(ga::normalize(b)), "normalize(line)");
        check(any_nonfinite(ga::inverse(b)), "inverse(line)");
        check(any_nonfinite(ga::dual(b)), "dual(line)");
        check(any_nonfinite(ga::exp(b)), "exp(line)");
        check(any_nonfinite(b * l2), "line * line");
        check(any_nonfinite(ga::meet(w, b)), "meet(plane, NaN line)");
        check(any_nonfinite(ga::join(q, b)), "join(point, NaN line)");
        check(any_nonfinite(ga::sandwich(b, m2)), "sandwich(line, motor)");
        check(any_nonfinite(ga::project(b, w)), "project(line, plane)");
        check(any_nonfinite(ga::project(q, b)), "project(point, NaN line)");
        check(any_nonfinite(ga::distance(q, b)), "distance(point, NaN line)");
        check(any_nonfinite(ga::motor_from_line_angle(b, T(0.5))), "motor_from_line_angle(NaN line)");
        check(any_nonfinite(ga::screw(b, T(0.5), T(1))), "screw(NaN line)");
    }

    // Points: e012..e023 ideal, e123 the weight.
    for (int i = 0; i < 4; ++i) {
        Trivector<T> p = ga::point(T(1), T(2), T(3));
        T* const slot[4] = {&p.e012, &p.e013, &p.e023, &p.e123};
        *slot[i] = nan;
        if (i == 3) {
            check(any_nonfinite(ga::norm(p)), "norm(point)");
        }
        check(any_nonfinite(ga::normalize(p)), "normalize(point)");
        check(any_nonfinite(ga::inverse(p)), "inverse(point)");
        check(any_nonfinite(ga::dual(p)), "dual(point)");
        check(any_nonfinite(ga::position(p)), "position(point)");
        check(any_nonfinite(ga::join(p, q)), "join(point, point)");
        check(any_nonfinite(ga::join(p, l2)), "join(point, line)");
        check(any_nonfinite(ga::sandwich(p, m2)), "sandwich(point, motor)");
        check(any_nonfinite(ga::project(p, w)), "project(point, plane)");
        check(any_nonfinite(ga::project(p, l2)), "project(point, line)");
        check(any_nonfinite(ga::reflect(p, w)), "reflect(point, plane)");
        check(any_nonfinite(ga::distance(p, q)), "distance(point, point)");
        check(any_nonfinite(ga::distance(p, w)), "distance(point, plane)");
        check(any_nonfinite(ga::distance(p, l2)), "distance(point, line)");
        check(any_nonfinite(ga::translator_between(p, q)), "translator_between(NaN point, point)");
    }

    // Motors: e01..e03 and e0123 ideal, the rest Euclidean.
    for (int i = 0; i < 8; ++i) {
        Motor<T> m = m2;
        T* const slot[8] = {&m.scalar.value, &m.bivector.e01, &m.bivector.e02, &m.bivector.e03,
                            &m.bivector.e12, &m.bivector.e13, &m.bivector.e23, &m.quadvector.e0123};
        *slot[i] = nan;
        const bool euclidean = i == 0 || (i >= 4 && i <= 6);
        if (euclidean) {
            check(any_nonfinite(ga::norm(m)), "norm(motor)");
        }
        check(any_nonfinite(ga::normalize(m)), "normalize(motor)");
        check(any_nonfinite(ga::inverse(m)), "inverse(motor)");
        check(any_nonfinite(ga::reverse(m)), "reverse(motor)");
        check(any_nonfinite(ga::log(m)), "log(motor)");
        check(any_nonfinite(m * m2), "motor * motor");
        check(any_nonfinite(ga::sandwich(q, m)), "sandwich(point, NaN motor)");
        check(any_nonfinite(ga::sandwich(l2, m)), "sandwich(line, NaN motor)");
        check(any_nonfinite(ga::slerp(m, m2, T(0.5))), "slerp(NaN, m, 0.5)");
        check(any_nonfinite(ga::slerp(m2, m, T(0.5))), "slerp(m, NaN, 0.5)");
    }

    // A NaN parameter rather than a NaN operand.
    check(any_nonfinite(ga::rotor_from_axis_angle(T(1), T(0), T(0), nan)), "rotor_from_axis_angle(axis, NaN)");
    check(any_nonfinite(ga::rotor_from_axis_angle(nan, T(0), T(0), T(1))), "rotor_from_axis_angle(NaN axis)");
    check(any_nonfinite(ga::motor_from_line_angle(l2, nan)), "motor_from_line_angle(line, NaN)");
    check(any_nonfinite(ga::screw(l2, T(0.5), nan)), "screw(line, angle, NaN)");
    check(any_nonfinite(ga::translator(nan, T(0), T(0))), "translator(NaN)");
    check(any_nonfinite(ga::point(nan, T(0), T(0))), "point(NaN)");
    check(any_nonfinite(ga::slerp(ga::identity_motor<T>(), m2, nan)), "slerp(a, b, NaN)");
    check(any_nonfinite(ga::exp(Bivector<T>(T(0), T(0), T(0), std::numeric_limits<T>::infinity(), T(0), T(0)))),
          "exp(infinite line)");
}

// norm() ignores the ideal components by design, and that includes a NaN in
// one: a poisoned plane can report a perfectly ordinary norm. Pinned so the
// blind spot is documented by a test, not just by docs/safety-audit.md (A4).
void test_norm_ignores_ideal_nan()
{
    section("norm ignores a NaN in an ideal component (docs/safety-audit.md A4)");

    const double nan = std::numeric_limits<double>::quiet_NaN();
    check(ga::norm(Vector<double>(nan, 1.0, 0.0, 0.0)).value == 1.0, "norm(plane with NaN e0) is 1");
    check(ga::norm(Bivector<double>(nan, 0.0, 0.0, 1.0, 0.0, 0.0)).value == 1.0, "norm(line with NaN e01) is 1");
    check(ga::norm(Trivector<double>(nan, 0.0, 0.0, 1.0)).value == 1.0, "norm(point with NaN e012) is 1");
    check((Vector<double>(nan, 1.0, 0.0, 0.0) | Vector<double>(0.0, 1.0, 0.0, 0.0)).value == 1.0,
          "plane | plane with NaN e0 is finite");
    check(!ga::is_ideal(Trivector<double>(0.0, 0.0, 0.0, nan)), "is_ideal(point with NaN weight) is false");
}

// ---------------------------------------------------------------------------
// Zero and ideal inputs: today's behaviour, pinned
// ---------------------------------------------------------------------------

void test_zero_inputs()
{
    section("zero and ideal inputs (today's behaviour; see docs/safety-audit.md A1, A3)");

    const Vector<double> z1 = ga::plane(0.0, 0.0, 1.0, -1.0);
    const Trivector<double> idealPoint = ga::ideal_point(1.0, 2.0, 3.0);
    const Bivector<double> idealLine(1.0, 2.0, 3.0, 0.0, 0.0, 0.0);

    // A1: a degenerate input returns a value that looks like a real answer.
    check(ga::distance(idealPoint, z1) == 0.0, "distance(ideal point, plane) is 0, as if on the plane");
    check(ga::distance(ga::point(0.0, 0.0, 5.0), ga::ideal_plane<double>()) == 0.0,
          "distance(point, ideal plane) is 0");
    check(ga::distance(idealPoint, ga::point(1.0, 2.0, 3.0)) == 0.0, "distance(ideal point, point) is 0");
    check(ga::distance(ga::point(1.0, 2.0, 3.0), idealLine) == 0.0, "distance(point, ideal line) is 0");
    const std::array<double, 3> pos = ga::position(idealPoint);
    check(pos[0] == 0.0 && pos[1] == 0.0 && pos[2] == 0.0, "position(ideal point) is the origin");
    check(ga::project(idealLine, ga::ideal_plane<double>()) == Bivector<double>(),
          "project(line, ideal plane) is the zero line");
    check(ga::log(Motor<double>()) == Bivector<double>(), "log(zero motor) is the zero line");
    check(ga::rotor_from_axis_angle(0.0, 0.0, 0.0, 1.0) == ga::identity_rotor<double>(),
          "rotor_from_axis_angle(zero axis) is the identity rotor");
    check(ga::motor_between(ga::ideal_plane<double>(), z1) == ga::identity_motor<double>(),
          "motor_between(ideal plane, plane) is the identity motor");

    // A3: PGA's inverses return zero for a non-invertible input...
    check(ga::inverse(Vector<double>()) == Vector<double>(), "inverse(zero plane) is zero");
    check(ga::inverse(ga::ideal_plane<double>()) == Vector<double>(), "inverse(ideal plane) is zero");
    check(ga::inverse(idealLine) == Bivector<double>(), "inverse(ideal line) is zero");
    check(ga::inverse(idealPoint) == Trivector<double>(), "inverse(ideal point) is zero");
    check(ga::inverse(Multivector<double>()) == Multivector<double>(), "inverse(Multivector 0) is zero");
    check(ga::inverse(Rotor<double>(Scalar<double>(0.0), 0.0, 0.0, 0.0)) == Rotor<double>(),
          "inverse(zero rotor) is zero");
    check(ga::inverse(Translator<double>(Scalar<double>(0.0), 1.0, 2.0, 3.0)) == Translator<double>(),
          "inverse(ideal translator) is zero");
    // ...except the Scalar one, which divides. The zero goes through a
    // volatile so the optimizer cannot fold it into a constant division,
    // which MSVC at /O2 warns about (C4723) inside the library.
    volatile double opaque = 0.0;
    check(std::isinf(ga::inverse(Scalar<double>(opaque)).value), "inverse(Scalar 0) is infinite");
}

// ---------------------------------------------------------------------------
// The edges of the input domain
// ---------------------------------------------------------------------------

template<typename T>
void test_domain_edge(T s)
{
    // inverse(line) divides by the square of a squared norm: the tightest
    // expression in the library, and the one the domain was derived from.
    const Bivector<T> line(s, T(0), T(0), s, T(0), T(0));
    const Bivector<T> inv = ga::inverse(line);
    check_close_t(inv.e01 * s, T(-1), "inverse(line) ideal part at the edge");
    check_close_t(inv.e12 * s, T(-1), "inverse(line) Euclidean part at the edge");

    const Vector<T> unit = ga::normalize(ga::plane(s, T(2) * s, T(2) * s, s));
    check_close_t(unit.e1, T(1) / T(3), "normalize(plane) at the edge");
    check_close_t(unit.e0, T(1) / T(3), "normalize(plane) offset at the edge");

    check_close_t(ga::distance(ga::point(s, T(0), T(0)), ga::point(T(0), T(0), T(0))) / s, T(1),
                  "distance(point, point) at the edge");
    check_close_t(ga::distance(ga::point(T(0), T(0), s), ga::plane(T(0), T(0), T(1), T(0))) / s, T(1),
                  "distance(point, plane) at the edge");

    // A quarter turn about z carries (s, 0, 0) to (0, s, 0).
    const T halfPi = std::atan2(T(1), T(0));
    const std::array<T, 3> turned =
        ga::position(ga::sandwich(ga::point(s, T(0), T(0)), ga::rotor_from_axis_angle(T(0), T(0), T(1), halfPi)));
    check_close_t(turned[1] / s, T(1), "rotate a point at the edge: x goes to y");
    check_close_t(turned[0] / s + T(1), T(1), "rotate a point at the edge: no x left");

    const std::array<T, 3> moved =
        ga::position(ga::sandwich(ga::point(s, T(0), T(0)), ga::translator(s, T(0), T(0))));
    check_close_t(moved[0] / s, T(2), "translate a point at the edge");
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
    test_norm_ignores_ideal_nan();
    test_zero_inputs();
}

} // namespace cctest
