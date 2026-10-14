// Cl(3,0,1) tests -- properties: identities, round trips and invariants run in
// float, double and long double, with tolerances scaled to each type's epsilon.
//
// One category of the PGA suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

// ---------------------------------------------------------------------------
// Properties in every floating-point type
//
// Every other check in this suite runs in double, against values picked for
// double. These run the library in float, double and long double, and assert
// only what holds at any precision -- identities, round trips and invariants --
// with a tolerance scaled to the type's epsilon. float is where a precision bug
// shows first. (long double is 80-bit with MinGW, and the same as double with
// MSVC.)
// ---------------------------------------------------------------------------

// The sixteen components in the same order as make_mv and the docs.
template<typename T>
Multivector<T> make_mv_t(const T (&c)[16])
{
    return Multivector<T>(Scalar<T>(c[0]), Vector<T>(c[1], c[2], c[3], c[4]),
                          Bivector<T>(c[5], c[6], c[7], c[8], c[9], c[10]),
                          Trivector<T>(c[11], c[12], c[13], c[14]), Quadvector<T>(c[15]));
}

template<typename T>
Multivector<T> basis_blade_t(int k)
{
    T c[16] = {};
    c[k] = T(1);
    return make_mv_t(c);
}

template<typename T>
void check_vector_t(const Vector<T>& a, const Vector<T>& b, const std::string& what, T scale = T(1))
{
    check_close_t(a.e0, b.e0, what + ".e0", scale);
    check_close_t(a.e1, b.e1, what + ".e1", scale);
    check_close_t(a.e2, b.e2, what + ".e2", scale);
    check_close_t(a.e3, b.e3, what + ".e3", scale);
}

template<typename T>
void check_bivector_t(const Bivector<T>& a, const Bivector<T>& b, const std::string& what, T scale = T(1))
{
    check_close_t(a.e01, b.e01, what + ".e01", scale);
    check_close_t(a.e02, b.e02, what + ".e02", scale);
    check_close_t(a.e03, b.e03, what + ".e03", scale);
    check_close_t(a.e12, b.e12, what + ".e12", scale);
    check_close_t(a.e13, b.e13, what + ".e13", scale);
    check_close_t(a.e23, b.e23, what + ".e23", scale);
}

template<typename T>
void check_trivector_t(const Trivector<T>& a, const Trivector<T>& b, const std::string& what, T scale = T(1))
{
    check_close_t(a.e012, b.e012, what + ".e012", scale);
    check_close_t(a.e013, b.e013, what + ".e013", scale);
    check_close_t(a.e023, b.e023, what + ".e023", scale);
    check_close_t(a.e123, b.e123, what + ".e123", scale);
}

template<typename T>
void check_motor_t(const Motor<T>& a, const Motor<T>& b, const std::string& what, T scale = T(1))
{
    check_close_t(a.scalar.value, b.scalar.value, what + ".scalar", scale);
    check_bivector_t(a.bivector, b.bivector, what + ".bivector", scale);
    check_close_t(a.quadvector.e0123, b.quadvector.e0123, what + ".e0123", scale);
}

template<typename T>
void check_mv_t(const Multivector<T>& a, const Multivector<T>& b, const std::string& what, T scale = T(1))
{
    check_close_t(a.scalar.value, b.scalar.value, what + ".scalar", scale);
    check_vector_t(a.vector, b.vector, what + ".vector", scale);
    check_bivector_t(a.bivector, b.bivector, what + ".bivector", scale);
    check_trivector_t(a.trivector, b.trivector, what + ".trivector", scale);
    check_close_t(a.quadvector.e0123, b.quadvector.e0123, what + ".e0123", scale);
}

template<typename T>
void test_properties()
{
    static const std::string label = std::string("properties<") + type_name<T>() + ">";
    section(label.c_str());

    const T pi = T(3.141592653589793238462643383279502884L);
    const Multivector<T> one = basis_blade_t<T>(0);
    const Motor<T> identity = ga::identity_motor<T>();

    // The table is exact in every type: each basis coefficient is 0 or +-1,
    // and the 64 products containing e0 twice are exactly zero.
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            int k = -1;
            const int sign = reference_blade_product(i, j, k);
            const Multivector<T> expected = sign == 0 ? Multivector<T>()
                                                      : basis_blade_t<T>(k) * Scalar<T>(T(sign));
            check(basis_blade_t<T>(i) * basis_blade_t<T>(j) == expected,
                  std::string("cayley ") + kBladeNames[i] + " * " + kBladeNames[j]);
        }
    }

    // Algebraic laws, on multivectors with few-bit components.
    const T ca[16] = {0.5, -1.25, 2, 0.75, 1, -0.5, 0.25, -1.5, 1.25, 0.5, -2, 0.75, -0.25, 1, 2, -0.5};
    const T cb[16] = {-1.5, 0.25, -2, 1, 0.5, 2, -1, 0.75, -0.5, 1.5, 0.25, -1, 2, -0.75, 0.5, 1};
    const T cc[16] = {2, 0.5, -0.75, -1.25, -2, 1, 0.5, 0.25, 1, -1.5, 0.75, 2, -0.5, 0.25, -1, 1.5};
    const Multivector<T> a = make_mv_t(ca), b = make_mv_t(cb), c = make_mv_t(cc);
    check_mv_t((a * b) * c, a * (b * c), "associativity");
    check_mv_t(a * (b + c), a * b + a * c, "distributivity");
    check_mv_t(ga::reverse(a * b), ga::reverse(b) * ga::reverse(a), "reverse(ab) == ~b ~a");

    // Every invertible type. The line and the motor both carry an e0 part, so
    // the study-number correction in their inverses is exercised.
    const Vector<T> p = ga::plane<T>(1, -2, 2, 3);
    const Bivector<T> line = ga::line_through_points(ga::point<T>(1, 2, 3), ga::point<T>(-1, 0, 2));
    const Trivector<T> pt = ga::point<T>(1, -2, 3) * Scalar<T>(T(2.5));
    const Rotor<T> r = ga::rotor_from_axis_angle<T>(1, 2, -2, T(0.7)) * Scalar<T>(3);
    const Translator<T> tr = ga::translator<T>(1, -2, 3) * Scalar<T>(2);
    const Motor<T> screw = ga::screw(ga::normalize(line), T(0.9), T(1.5));
    const Motor<T> drifted = Motor<T>(screw.scalar, screw.bivector,
                                      screw.quadvector + Quadvector<T>(T(0.01))) * Scalar<T>(T(1.3));
    check_mv_t(p * ga::inverse(p), one, "plane * inverse(plane)");
    check_mv_t(line * ga::inverse(line), one, "line * inverse(line)");
    check_mv_t(pt * ga::inverse(pt), one, "point * inverse(point)");
    check_mv_t(ga::to_multivector(r * ga::inverse(r)), one, "rotor * inverse(rotor)");
    check_mv_t(ga::to_multivector(tr * ga::inverse(tr)), one, "translator * inverse(translator)");
    check_motor_t(drifted * ga::inverse(drifted), identity, "drifted motor * inverse");
    check_motor_t(ga::inverse(drifted) * drifted, identity, "inverse * drifted motor");
    const Motor<T> renormalized = ga::normalize(drifted);
    check_motor_t(renormalized * ga::reverse(renormalized), identity, "normalize(drifted) is unit, e0123 included");
    const Multivector<T> m = a;
    check_mv_t(m * ga::inverse(m), one, "multivector * inverse", T(4));
    check_mv_t(ga::inverse(m) * m, one, "inverse * multivector", T(4));

    // exp and log of a screw -- a rotation with a slide along its axis -- across
    // the range and close to the identity.
    const T angles[] = {T(-1.5), T(-0.25), T(-1e-2), T(-1e-3), T(1e-3), T(1e-2), T(0.25), T(1.5)};
    for (const T theta : angles) {
        const std::string at = " at " + describe(theta);
        // Euclidean part theta * a unit direction; an ideal part that is not a
        // multiple of it, so the result is a genuine screw with an offset axis.
        const Bivector<T> bt(T(0.3), T(-0.2), T(0.5), T(0.6) * theta, T(0) * theta, T(0.8) * theta);
        const Motor<T> e = ga::exp(bt);
        check_motor_t(e * ga::reverse(e), identity, "exp is a unit motor" + at);
        check_bivector_t(ga::log(e), bt, "log(exp(b)) == b" + at);
        check_motor_t(ga::exp(ga::log(e)), e, "exp(log(m)) == m" + at);
    }

    // Rigid motions: exact translation, the right-hand rule, a half turn about
    // an offset axis, and distances preserved.
    check_trivector_t(ga::sandwich(ga::point<T>(1, 2, 3), ga::translator<T>(4, -1, 2)), ga::point<T>(5, 1, 5),
                      "translator moves a point");
    check_trivector_t(ga::sandwich(ga::point<T>(1, 0, 0), ga::rotor_from_axis_angle<T>(0, 0, 1, pi / T(2))),
                      ga::point<T>(0, 1, 0), "+90 about +z takes (1,0,0) to (0,1,0)");
    const Bivector<T> offset_axis = ga::line_through_points(ga::point<T>(1, 0, 0), ga::point<T>(1, 0, 1));
    check_trivector_t(ga::sandwich(ga::point<T>(2, 0, 0), ga::motor_from_line_angle(offset_axis, pi)),
                      ga::point<T>(0, 0, 0), "half turn about the line x = 1 takes (2,0,0) to the origin");
    const Trivector<T> P = ga::point<T>(1, -2, 3), Q = ga::point<T>(-0.5, 4, 1);
    check_close_t(ga::distance(ga::sandwich(P, screw), ga::sandwich(Q, screw)), ga::distance(P, Q),
                  "a screw preserves distances");
    check_trivector_t(ga::sandwich(P, ga::translator_between(P, Q)), Q, "translator_between lands");

    // motor_between lands, including on nearly opposite planes.
    const Vector<T> from = ga::plane<T>(1, 2, -2, 3);
    const Bivector<T> hinge = ga::line_through_origin<T>(2, -1, 0);   // perpendicular to from's normal
    const Vector<T> targets[] = {
        ga::plane<T>(0, 1, 0, -2), ga::plane<T>(-1, 1, 1, 0.5),
        ga::sandwich(-from, ga::motor_from_line_angle(hinge, T(0.5))),
        ga::sandwich(-from, ga::motor_from_line_angle(hinge, T(1e-2))),
        ga::sandwich(-from, ga::motor_from_line_angle(hinge, T(1e-3))),
        ga::sandwich(-from, ga::motor_from_line_angle(hinge, T(1e-5))),
        ga::plane<T>(-1, -2, 2, 5),   // exactly opposite normal, different offset
    };
    for (const Vector<T>& to : targets) {
        const std::string at = " to (" + describe(to.e0) + ", " + describe(to.e1) + ", " + describe(to.e2)
                             + ", " + describe(to.e3) + ")";
        const Motor<T> between = ga::motor_between(from, to);
        check_motor_t(between * ga::reverse(between), identity, "motor_between is unit" + at);
        check_vector_t(ga::sandwich(ga::normalize(from), between), ga::normalize(to), "motor_between lands" + at);
    }

    check_motor_t(ga::slerp(identity, screw, T(0)), identity, "slerp at 0");
    check_motor_t(ga::slerp(identity, screw, T(1)), screw, "slerp at 1");

    // Angles between planes and between lines, including nearly parallel ones.
    const Bivector<T> hinge2 = ga::line_through_origin<T>(1, -1, 0);   // perpendicular to line's direction
    const T alphas[] = {T(1e-3), T(1e-2), T(0.3), T(1.2), T(2.8)};
    for (const T alpha : alphas) {
        const Vector<T> turned = ga::sandwich(from, ga::motor_from_line_angle(hinge, alpha));
        check_close_t(ga::angle(from, turned), alpha, "angle between planes at " + describe(alpha));
        const Bivector<T> turned_line = ga::sandwich(line, ga::motor_from_line_angle(hinge2, alpha));
        check_close_t(ga::angle(line, turned_line), alpha, "angle between lines at " + describe(alpha));
    }

    // The default tolerance of approx_equal absorbs a few ulps, not more.
    const T eps = std::numeric_limits<T>::epsilon();
    check(ga::approx_equal(Vector<T>(1, 2, 3, 4), Vector<T>(T(1) + T(4) * eps, 2, 3, 4)),
          "approx_equal accepts 4 ulps");
    check(!ga::approx_equal(Vector<T>(1, 2, 3, 4), Vector<T>(T(1) + T(1000) * eps, 2, 3, 4)),
          "approx_equal rejects 1000 ulps");
}

void test_properties_in_every_type()
{
    test_properties<float>();
    test_properties<double>();
    test_properties<long double>();
}

} // namespace

void run_properties()
{
    test_properties_in_every_type();
}

} // namespace cctest
