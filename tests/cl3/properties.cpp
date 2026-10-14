// Cl(3,0) tests -- properties: identities, round trips and invariants run in
// float, double and long double, with tolerances scaled to each type's epsilon.
//
// One category of the Cl3 suite; main.cpp runs them all, and support.hpp
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
Multivector<T> make_mv_t(T s, T x, T y, T z, T xy, T xz, T yz, T e123)
{
    return Multivector<T>(Scalar<T>(s), Vector<T>(x, y, z), Bivector<T>(xy, xz, yz), Trivector<T>(e123));
}

// The eight basis blades in the library's order, as index lists, and an
// independent reference for their products: concatenate the lists, bubble-sort
// counting swaps, then cancel equal neighbours (e_i^2 = +1 never flips a sign).
struct BladeIndices
{
    int count;
    int index[3];
};

const BladeIndices kBlades[8] = {
    {0, {0, 0, 0}}, {1, {1, 0, 0}}, {1, {2, 0, 0}}, {1, {3, 0, 0}},
    {2, {1, 2, 0}}, {2, {1, 3, 0}}, {2, {2, 3, 0}}, {3, {1, 2, 3}},
};

const char* const kBladeNames[8] = {"1", "e1", "e2", "e3", "e12", "e13", "e23", "e123"};

int reference_blade_product(int i, int j, int& result_blade)
{
    int list[6];
    int n = 0;
    for (int k = 0; k < kBlades[i].count; ++k) list[n++] = kBlades[i].index[k];
    for (int k = 0; k < kBlades[j].count; ++k) list[n++] = kBlades[j].index[k];

    int sign = 1;
    bool changed = true;
    while (changed) {
        changed = false;
        for (int k = 0; k + 1 < n; ++k) {
            if (list[k] > list[k + 1]) {
                const int tmp = list[k]; list[k] = list[k + 1]; list[k + 1] = tmp;
                sign = -sign;
                changed = true;
            }
        }
    }

    int out[3];
    int m = 0;
    for (int k = 0; k < n;) {
        if (k + 1 < n && list[k] == list[k + 1]) {
            k += 2;
        } else {
            out[m++] = list[k];
            k += 1;
        }
    }

    result_blade = -1;
    for (int b = 0; b < 8; ++b) {
        if (kBlades[b].count != m) continue;
        bool same = true;
        for (int k = 0; k < m; ++k) same = same && kBlades[b].index[k] == out[k];
        if (same) { result_blade = b; break; }
    }
    return sign;
}

template<typename T>
Multivector<T> basis_blade_t(int k)
{
    T c[8] = {T(0), T(0), T(0), T(0), T(0), T(0), T(0), T(0)};
    c[k] = T(1);
    return make_mv_t(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]);
}

template<typename T>
void check_vector_t(const Vector<T>& a, const Vector<T>& b, const std::string& what, T scale = T(1))
{
    check_close_t(a.x, b.x, what + ".x", scale);
    check_close_t(a.y, b.y, what + ".y", scale);
    check_close_t(a.z, b.z, what + ".z", scale);
}

template<typename T>
void check_bivector_t(const Bivector<T>& a, const Bivector<T>& b, const std::string& what, T scale = T(1))
{
    check_close_t(a.xy, b.xy, what + ".xy", scale);
    check_close_t(a.xz, b.xz, what + ".xz", scale);
    check_close_t(a.yz, b.yz, what + ".yz", scale);
}

template<typename T>
void check_rotor_t(const Rotor<T>& a, const Rotor<T>& b, const std::string& what, T scale = T(1))
{
    check_close_t(a.scalar.value, b.scalar.value, what + ".scalar", scale);
    check_bivector_t(a.bivector, b.bivector, what + ".bivector", scale);
}

template<typename T>
void check_mv_t(const Multivector<T>& a, const Multivector<T>& b, const std::string& what, T scale = T(1))
{
    check_close_t(a.scalar.value, b.scalar.value, what + ".scalar", scale);
    check_vector_t(a.vector, b.vector, what + ".vector", scale);
    check_bivector_t(a.bivector, b.bivector, what + ".bivector", scale);
    check_close_t(a.trivector.e123, b.trivector.e123, what + ".e123", scale);
}

template<typename T>
void test_properties()
{
    static const std::string label = std::string("properties<") + type_name<T>() + ">";
    section(label.c_str());

    const T pi = T(3.141592653589793238462643383279502884L);
    const Multivector<T> one = make_mv_t<T>(1, 0, 0, 0, 0, 0, 0, 0);
    const Rotor<T> identity = ga::identity_rotor<T>();

    // The table is exact in every type: each basis coefficient is 0 or +-1.
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 8; ++j) {
            int k = -1;
            const int sign = reference_blade_product(i, j, k);
            check(basis_blade_t<T>(i) * basis_blade_t<T>(j)
                      == basis_blade_t<T>(k) * Scalar<T>(T(sign)),
                  std::string("cayley ") + kBladeNames[i] + " * " + kBladeNames[j]);
        }
    }

    // Algebraic laws, on multivectors with few-bit components.
    const Multivector<T> a = make_mv_t<T>(0.5, -1.25, 2, 0.75, 1, -0.5, 0.25, -1.5);
    const Multivector<T> b = make_mv_t<T>(-1.5, 0.25, -2, 1, 0.5, 2, -1, 0.75);
    const Multivector<T> c = make_mv_t<T>(2, 0.5, -0.75, -1.25, -2, 1, 0.5, 0.25);
    check_mv_t((a * b) * c, a * (b * c), "associativity");
    check_mv_t(a * (b + c), a * b + a * c, "distributivity");
    check_mv_t(ga::reverse(a * b), ga::reverse(b) * ga::reverse(a), "reverse(ab) == ~b ~a");
    check_mv_t(ga::conjugate(a * b), ga::conjugate(b) * ga::conjugate(a), "conjugate(ab)");

    // Every invertible type.
    const Vector<T> v(3, -4, 12);
    const Bivector<T> bv(T(2.5), -1, 2);
    const Trivector<T> tv(T(-1.5));
    const Rotor<T> r(Scalar<T>(3), Bivector<T>(4, -2, 1));
    const Multivector<T> m = make_mv_t<T>(1, -2, 3, 0.5, -1, 2, 0.25, 1.5);
    check_mv_t(v * ga::inverse(v), one, "v * inverse(v)");
    check_mv_t(bv * ga::inverse(bv), one, "b * inverse(b)");
    check_close_t((tv * ga::inverse(tv)).value, T(1), "t * inverse(t)");
    check_rotor_t(r * ga::inverse(r), identity, "r * inverse(r)");
    check_mv_t(m * ga::inverse(m), one, "m * inverse(m)");
    check_mv_t(ga::inverse(m) * m, one, "inverse(m) * m");

    check_close_t(ga::norm(ga::normalize(v)).value, T(1), "normalize(Vector) is unit");
    check_close_t(ga::norm(ga::normalize(bv)).value, T(1), "normalize(Bivector) is unit");
    check_close_t(ga::norm(ga::normalize(r)).value, T(1), "normalize(Rotor) is unit");
    check_close_t(ga::norm(ga::normalize(m)).value, T(1), "normalize(Multivector) is unit");

    // exp and log, across the open range and close to the identity, in a
    // plane that is not a coordinate plane.
    const Bivector<T> plane = ga::normalize(Bivector<T>(1, 2, -2));
    const T angles[] = {T(-3), T(-1.5), T(-0.25), T(-1e-2), T(-1e-3), T(-1e-4),
                        T(1e-4), T(1e-3), T(1e-2), T(0.25), T(1.5), T(3)};
    for (const T theta : angles) {
        const std::string at = " at " + describe(theta);
        const Bivector<T> bt = plane * Scalar<T>(theta);
        const Rotor<T> e = ga::exp(bt);
        check_close_t(ga::norm(e).value, T(1), "exp is unit" + at);
        check_bivector_t(ga::log(e), bt, "log(exp(b)) == b" + at);
        check_rotor_t(ga::exp(ga::log(e)), e, "exp(log(r)) == r" + at);
    }

    // Right-hand rule, and rotations preserve lengths and angles.
    const Vector<T> e1(1, 0, 0), e2(0, 1, 0), e3(0, 0, 1);
    check_vector_t(ga::rotate(e1, ga::rotor_from_axis_angle(e3, pi / T(2))), e2, "+90 about +z takes e1 to e2");
    const Vector<T> w(T(-0.5), 3, T(1.25));
    const Vector<T> axis = ga::normalize(Vector<T>(1, -1, 2));
    const Rotor<T> ru = ga::normalize(r);
    for (int i = -31; i <= 31; ++i) {
        const Rotor<T> q = ga::rotor_from_axis_angle(axis, T(0.1) * T(i));
        const std::string at = " at " + describe(T(0.1) * T(i));
        check_close_t(ga::norm(ga::rotate(v, q)).value, ga::norm(v).value, "rotation keeps length" + at);
        check_close_t((ga::rotate(v, q) | ga::rotate(w, q)).value, (v | w).value, "rotation keeps dot" + at);
        check_vector_t(ga::rotate(ga::rotate(v, q), ru), ga::rotate(v, ru * q), "composition" + at, T(4));
    }

    // rotor_between lands, including on nearly opposite vectors -- the case
    // where 1 + a.b cancels catastrophically.
    const Vector<T> from(3, 1, -2);
    const Vector<T> across = ga::normalize(Vector<T>(1, -1, 1));   // perpendicular to `from`
    const Vector<T> targets[] = {
        Vector<T>(-1, 2, 0), Vector<T>(0, -5, 1), Vector<T>(T(2.5), T(0.1), T(-2)),
        ga::rotate(-from, ga::rotor_from_axis_angle(across, T(0.5))),
        ga::rotate(-from, ga::rotor_from_axis_angle(across, T(1e-2))),
        ga::rotate(-from, ga::rotor_from_axis_angle(across, T(1e-3))),
        ga::rotate(-from, ga::rotor_from_axis_angle(across, T(1e-5))),
    };
    for (const Vector<T>& to : targets) {
        const std::string at = " to (" + describe(to.x) + ", " + describe(to.y) + ", " + describe(to.z) + ")";
        const Rotor<T> between = ga::rotor_between(from, to);
        check_close_t(ga::norm(between).value, T(1), "rotor_between is unit" + at);
        check_vector_t(ga::rotate(ga::normalize(from), between), ga::normalize(to),
                       "rotor_between lands" + at);
    }
    const Rotor<T> half = ga::rotor_between(from, -from);
    check_vector_t(ga::rotate(ga::normalize(from), half), -ga::normalize(from), "rotor_between of opposites");

    const Rotor<T> ra = ga::rotor_from_axis_angle(axis, T(0.2));
    const Rotor<T> rb = ga::rotor_from_axis_angle(axis, T(1.0));
    check_rotor_t(ga::slerp(ra, rb, T(0)), ra, "slerp at 0");
    check_rotor_t(ga::slerp(ra, rb, T(1)), rb, "slerp at 1");
    check_rotor_t(ga::slerp(ra, rb, T(0.5)), ga::rotor_from_axis_angle(axis, (T(0.2) + T(1.0)) / T(2)),
                  "slerp midpoint");

    // Geometry identities.
    const Vector<T> n(1, 2, -2);
    check_vector_t(ga::reflect(ga::reflect(v, n), n), v, "reflecting twice is the identity");
    check_vector_t(ga::project(v, n) + ga::reject(v, n), v, "project + reject == v");
    check_vector_t(ga::project(v, bv) + ga::reject(v, bv), v, "project + reject == v, plane");

    // The default tolerance of approx_equal absorbs a few ulps, not more.
    const T eps = std::numeric_limits<T>::epsilon();
    check(ga::approx_equal(Vector<T>(1, 2, 3), Vector<T>(T(1) + T(4) * eps, 2, 3)),
          "approx_equal accepts 4 ulps");
    check(!ga::approx_equal(Vector<T>(1, 2, 3), Vector<T>(T(1) + T(1000) * eps, 2, 3)),
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
