// Cl(2,0) tests -- properties: identities, round trips and invariants run in
// float, double and long double, with tolerances scaled to each type's epsilon.
//
// One category of the Cl2 suite; main.cpp runs them all, and support.hpp
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

template<typename T>
Multivector<T> make_mv_t(T s, T x, T y, T xy)
{
    return Multivector<T>(Scalar<T>(s), Vector<T>(x, y), Bivector<T>(xy));
}

template<typename T>
Multivector<T> basis_blade_t(int k)
{
    T c[4] = {T(0), T(0), T(0), T(0)};
    c[k] = T(1);
    return make_mv_t(c[0], c[1], c[2], c[3]);
}

template<typename T>
void check_vector_t(const Vector<T>& a, const Vector<T>& b, const std::string& what, T scale = T(1))
{
    check_close_t(a.x, b.x, what + ".x", scale);
    check_close_t(a.y, b.y, what + ".y", scale);
}

template<typename T>
void check_rotor_t(const Rotor<T>& a, const Rotor<T>& b, const std::string& what, T scale = T(1))
{
    check_close_t(a.scalar.value, b.scalar.value, what + ".scalar", scale);
    check_close_t(a.bivector.xy, b.bivector.xy, what + ".bivector.xy", scale);
}

template<typename T>
void check_mv_t(const Multivector<T>& a, const Multivector<T>& b, const std::string& what, T scale = T(1))
{
    check_close_t(a.scalar.value, b.scalar.value, what + ".scalar", scale);
    check_vector_t(a.vector, b.vector, what + ".vector", scale);
    check_close_t(a.bivector.xy, b.bivector.xy, what + ".bivector.xy", scale);
}

template<typename T>
void test_properties()
{
    static const std::string label = std::string("properties<") + type_name<T>() + ">";
    section(label.c_str());

    const T pi = T(3.141592653589793238462643383279502884L);
    const Multivector<T> one = make_mv_t<T>(1, 0, 0, 0);
    const Rotor<T> identity = ga::identity_rotor<T>();

    // The table is exact in every type: each basis coefficient is 0 or +-1.
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            int k = -1;
            const int sign = reference_blade_product(i, j, k);
            const std::string what = std::string("cayley ") + kBladeNames[i] + " * " + kBladeNames[j];
            check(basis_blade_t<T>(i) * basis_blade_t<T>(j)
                      == basis_blade_t<T>(k) * Scalar<T>(T(sign)),
                  what);
        }
    }

    // Algebraic laws, on multivectors with few-bit components.
    const Multivector<T> a = make_mv_t<T>(0.5, -1.25, 2, 0.75);
    const Multivector<T> b = make_mv_t<T>(-1.5, 0.25, -2, 1);
    const Multivector<T> c = make_mv_t<T>(2, 0.5, -0.75, -1.25);
    check_mv_t((a * b) * c, a * (b * c), "associativity");
    check_mv_t(a * (b + c), a * b + a * c, "distributivity");
    check_mv_t(ga::reverse(a * b), ga::reverse(b) * ga::reverse(a), "reverse(ab) == ~b ~a");
    check_mv_t(ga::conjugate(a * b), ga::conjugate(b) * ga::conjugate(a), "conjugate(ab)");

    // Every invertible type, from both sides where the product is not closed.
    const Vector<T> v(3, -4);
    const Bivector<T> bv(T(2.5));
    const Rotor<T> r(Scalar<T>(3), Bivector<T>(4));
    const Multivector<T> m = make_mv_t<T>(1, -2, 3, 0.5);
    check_mv_t(v * ga::inverse(v), one, "v * inverse(v)");
    check_close_t((bv * ga::inverse(bv)).value, T(1), "b * inverse(b)");
    check_rotor_t(r * ga::inverse(r), identity, "r * inverse(r)");
    check_mv_t(m * ga::inverse(m), one, "m * inverse(m)");
    check_mv_t(ga::inverse(m) * m, one, "inverse(m) * m");

    check_close_t(ga::norm(ga::normalize(v)).value, T(1), "normalize(Vector) is unit");
    check_close_t(ga::norm(ga::normalize(r)).value, T(1), "normalize(Rotor) is unit");
    check_close_t(ga::norm(ga::normalize(m)).value, T(1), "normalize(Multivector) is unit");

    // exp and log, and the rotor angle, across the whole open range.
    for (int i = -31; i <= 31; ++i) {
        const T theta = T(0.1) * T(i);
        const std::string at = " at " + describe(theta);
        const Rotor<T> e = ga::exp(Bivector<T>(theta));
        check_close_t(ga::norm(e).value, T(1), "exp is unit" + at);
        check_close_t(ga::log(e).xy, theta, "log(exp(b)) == b" + at);
        check_rotor_t(ga::exp(ga::log(e)), e, "exp(log(r)) == r" + at);
        check_close_t(ga::rotor_angle(ga::rotor_from_angle(theta)), theta, "rotor_angle round trip" + at);
    }

    // Rotation preserves length and angles, and the sandwich is one-sided,
    // for rotors of any magnitude.
    const Vector<T> w(T(-0.5), T(3));
    for (int i = -31; i <= 31; ++i) {
        const Rotor<T> q = ga::rotor_from_angle(T(0.1) * T(i)) * Scalar<T>(T(1.5));
        const std::string at = " at " + describe(T(0.1) * T(i));
        const Rotor<T> u = ga::normalize(q);
        check_close_t(ga::norm(ga::rotate(v, u)).value, ga::norm(v).value, "rotation keeps length" + at);
        check_close_t(ga::angle_between(ga::rotate(v, u), ga::rotate(w, u)), ga::angle_between(v, w),
                      "rotation keeps angles" + at);
        check_vector_t(ga::sandwich(v, q), q * q * v, "R v ~R == R^2 v" + at);
    }

    // rotor_between lands, including on nearly opposite vectors -- the case
    // where 1 + a.b cancels catastrophically.
    const Vector<T> from(3, 1);
    const Vector<T> targets[] = {
        Vector<T>(-1, 2), Vector<T>(0, -5), Vector<T>(T(2.5), T(0.1)),
        ga::rotate(-from, ga::rotor_from_angle(T(0.5))),
        ga::rotate(-from, ga::rotor_from_angle(T(1e-1))),
        ga::rotate(-from, ga::rotor_from_angle(T(1e-2))),
        ga::rotate(-from, ga::rotor_from_angle(T(1e-3))),
        ga::rotate(-from, ga::rotor_from_angle(T(-1e-3))),
    };
    for (const Vector<T>& to : targets) {
        const std::string at = " to (" + describe(to.x) + ", " + describe(to.y) + ")";
        const Rotor<T> between = ga::rotor_between(from, to);
        check_close_t(ga::norm(between).value, T(1), "rotor_between is unit" + at);
        check_vector_t(ga::rotate(ga::normalize(from), between), ga::normalize(to),
                       "rotor_between lands" + at);
    }

    const Rotor<T> ra = ga::rotor_from_angle(T(0.2));
    const Rotor<T> rb = ga::rotor_from_angle(T(1.0));
    check_rotor_t(ga::slerp(ra, rb, T(0)), ra, "slerp at 0");
    check_rotor_t(ga::slerp(ra, rb, T(1)), rb, "slerp at 1");
    // Not T(0.6): the literals 0.2 and 1.0 are doubles, so in long double the
    // midpoint of the two angles is (T(0.2) + T(1.0)) / 2, which is not 0.6.
    check_close_t(ga::rotor_angle(ga::slerp(ra, rb, T(0.5))), (T(0.2) + T(1.0)) / T(2), "slerp midpoint");

    // Cl2 specifics: signed zero, perp, and the planar extras.
    check_close_t(ga::angle(-Vector<T>(1, 0)), pi, "angle(-e1) == +pi");
    check_close_t(ga::angle_between(Vector<T>(-1, 0), Vector<T>(1, 0)), pi, "angle_between(-e1, e1) == +pi");
    check_vector_t(ga::perp(v), ga::dual(v), "perp == dual");
    check_close_t(ga::signed_area(v, ga::perp(v)).value, ga::squared_norm(v).value,
                  "signed_area(v, perp v) == |v|^2");

    // The default tolerance of approx_equal absorbs a few ulps, not more.
    const T eps = std::numeric_limits<T>::epsilon();
    check(ga::approx_equal(Vector<T>(1, 2), Vector<T>(T(1) + T(4) * eps, 2)),
          "approx_equal accepts 4 ulps");
    check(!ga::approx_equal(Vector<T>(1, 2), Vector<T>(T(1) + T(1000) * eps, 2)),
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
