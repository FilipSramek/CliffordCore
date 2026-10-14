// Cl(2,0) tests -- rotations: exp and log, rotor construction, the sandwich
// product, and slerp.
//
// One category of the Cl2 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

// ---------------------------------------------------------------------------
// Rotations
// ---------------------------------------------------------------------------

void test_exp_log()
{
    section("exp and log");

    check_rotor(ga::exp(Bivector<double>(0.0)), 1.0, 0.0, "exp(0) is the identity");

    // Euler's formula, with e12 as i.
    for (int i = -30; i <= 30; ++i) {
        const double theta = 0.1 * i;
        const Rotor<double> r = ga::exp(Bivector<double>(theta));
        const std::string at = " at theta = " + std::to_string(theta);
        check_rotor(r, std::cos(theta), std::sin(theta), "exp(theta e12) == cos + sin e12" + at);
        check_close(ga::norm(r).value, 1.0, "exp is always unit" + at);
        check_bivector(ga::log(r), theta, "log(exp(b)) == b" + at);
    }

    // log undoes exp, and exp undoes log, on a unit rotor.
    const Rotor<double> r = ga::normalize(Rotor<double>(Scalar<double>(-0.3), Bivector<double>(0.9)));
    check_rotor_close(ga::exp(ga::log(r)), r, "exp(log(r)) == r");

    // A negative scalar part means a turn past a half: acos-based logs need a
    // clamp and a sign fix here, atan2 does not.
    check_bivector(ga::log(ga::exp(Bivector<double>(2.5))), 2.5, "log of a turn past a half");

    // atan2 only sees the ratio, so drift off unit length does not move the angle.
    const Rotor<double> drifted = r * Scalar<double>(3.7);
    check_bivector(ga::log(drifted), ga::log(r).xy, "log ignores a rotor's magnitude");
    check_bivector(ga::log(Rotor<double>()), 0.0, "log of the zero rotor is zero");

    // The range is (-pi, pi]: a -0.0 bivector part must not push -1 to -pi.
    check_bivector(ga::log(Rotor<double>(Scalar<double>(-1.0), Bivector<double>(-0.0))), kPi,
                   "log(-1 - 0.0 e12) == +pi e12");

    // exp(theta B) rotates by 2 theta -- clockwise for positive theta, which is
    // why rotor_from_angle carries a minus sign.
    const Vector<double> turned = ga::sandwich(Vector<double>(1, 0), ga::exp(Bivector<double>(0.4)));
    check_close(std::atan2(turned.y, turned.x), -0.8, "exp(0.4 e12) turns e1 by -0.8");
}

void test_rotor_construction()
{
    section("rotor construction");

    const Vector<double> e1(1, 0), e2(0, 1);

    check_rotor(ga::identity_rotor<double>(), 1.0, 0.0, "identity_rotor");
    check_vector_close(ga::rotate(Vector<double>(3, -2), ga::identity_rotor<double>()),
                       Vector<double>(3, -2), "the identity rotates nothing");

    // Counter-clockwise for a positive angle, the same sign as Cl(3,0) about +z.
    const Rotor<double> quarter = ga::rotor_from_angle(kPi / 2);
    check_rotor(quarter, std::cos(kPi / 4), -std::sin(kPi / 4), "rotor_from_angle(pi/2) components");
    check_vector_close(ga::rotate(e1, quarter), e2, "+90 degrees takes e1 to e2");
    check_vector_close(ga::rotate(e2, quarter), -e1, "+90 degrees takes e2 to -e1");
    check_vector_close(ga::rotate(e1, ga::rotor_from_angle(-kPi / 2)), -e2, "-90 degrees takes e1 to -e2");
    check_vector_close(ga::rotate(e1, ga::rotor_from_angle(kPi)), -e1, "180 degrees takes e1 to -e1");

    // Rotations add, and in two dimensions they commute.
    check_rotor_close(ga::rotor_from_angle(0.3) * ga::rotor_from_angle(0.9), ga::rotor_from_angle(1.2),
                      "rotor_from_angle(a) * rotor_from_angle(b) == rotor_from_angle(a + b)");

    // rotor_angle reads the angle back, in (-pi, pi].
    for (int i = -30; i <= 30; ++i) {
        const double theta = 0.1 * i;
        check_close(ga::rotor_angle(ga::rotor_from_angle(theta)), theta,
                    "rotor_angle round trip at " + std::to_string(theta));
    }
    const Rotor<double> r = ga::rotor_from_angle(1.1);
    check_close(ga::rotor_angle(r * Scalar<double>(-1.0)), 1.1, "rotor_angle(-r) == rotor_angle(r)");
    check_close(ga::rotor_angle(r * Scalar<double>(2.5)), 1.1, "rotor_angle ignores magnitude");
    check_close(ga::rotor_angle(ga::rotor_from_angle(1.5 * kPi)), -0.5 * kPi,
                "rotor_angle reports the equivalent angle in range");
    check_close(ga::rotor_angle(ga::identity_rotor<double>()), 0.0, "rotor_angle(identity) == 0");

    // Both exact half-turn rotors report +pi. For 0 + 1 e12 the numerator is
    // (-2 * 0) * 1 = -0.0, which plain atan2 would turn into -pi.
    check_close(ga::rotor_angle(Rotor<double>(Scalar<double>(0), Bivector<double>(1))), kPi,
                "rotor_angle(+e12) == +pi");
    check_close(ga::rotor_angle(Rotor<double>(Scalar<double>(0), Bivector<double>(-1))), kPi,
                "rotor_angle(-e12) == +pi");

    // rotor_between lands exactly, from non-unit inputs.
    const Vector<double> from(3.0, 1.0), to(-1.0, 2.0);
    const Rotor<double> between = ga::rotor_between(from, to);
    check_close(ga::norm(between).value, 1.0, "rotor_between is unit");
    check_vector_close(ga::rotate(ga::normalize(from), between), ga::normalize(to),
                       "rotor_between carries from onto to");
    check_close(ga::rotor_angle(ga::rotor_between(Vector<double>(1, 0), Vector<double>(0, 2))), kPi / 2,
                "rotor_between(e1, 2 e2) is a quarter turn");

    // Opposite: no plane to choose in 2D, only a direction. It is the exact
    // counter-clockwise half turn.
    const Rotor<double> half = ga::rotor_between(from, -from);
    check_rotor(half, 0.0, -1.0, "rotor_between of opposites is rotor_from_angle(pi)");
    check_vector_close(ga::rotate(from, half), -from, "and it lands on the opposite");

    check_rotor(ga::rotor_between(Vector<double>(), e1), 1.0, 0.0, "a zero input gives the identity");
}

void test_sandwich()
{
    section("sandwich");

    const Vector<double> v(2.0, -1.5);
    const double scales[3] = {0.5, 1.0, 2.5};

    for (int i = -30; i <= 30; ++i) {
        for (double k : scales) {
            // Non-unit rotors too: every identity below holds at any magnitude.
            const Rotor<double> r = ga::rotor_from_angle(0.1 * i) * Scalar<double>(k);
            const std::string at = " (angle " + std::to_string(0.1 * i) + ", scale "
                                 + std::to_string(k) + ")";

            // The closed form agrees with the two-sided general product.
            const Multivector<double> general =
                ga::to_multivector(r) * as_mv(v) * ga::to_multivector(ga::reverse(r));
            check_vector_close(ga::sandwich(v, r), general.vector, "closed form == R v ~R" + at, 1e-11);
            check_close(general.scalar.value, 0.0, "R v ~R has no scalar part" + at, 1e-11);
            check_close(general.bivector.xy, 0.0, "R v ~R has no bivector part" + at, 1e-11);

            // The 2D-only shortcut: the sandwich is a one-sided product.
            check_vector_close(ga::sandwich(v, r), r * r * v, "R v ~R == R^2 v" + at, 1e-11);
        }
    }

    // Rotation preserves length and the angle between two vectors.
    const Rotor<double> r = ga::rotor_from_angle(0.7);
    const Vector<double> w(-0.5, 3.0);
    check_close(ga::norm(ga::rotate(v, r)).value, ga::norm(v).value, "rotation preserves length");
    check_close(ga::angle_between(ga::rotate(v, r), ga::rotate(w, r)), ga::angle_between(v, w),
                "rotation preserves angles");
    check_close(ga::angle_between(v, ga::rotate(v, r)), 0.7, "rotate turns by the rotor's angle");

    // Composition: rotating by r1 then r2 is rotating by r2 * r1.
    const Rotor<double> r1 = ga::rotor_from_angle(0.4), r2 = ga::rotor_from_angle(-1.3);
    check_vector_close(ga::rotate(ga::rotate(v, r1), r2), ga::rotate(v, r2 * r1),
                       "rotate(rotate(v, r1), r2) == rotate(v, r2 r1)");

    // There is only one plane, so a unit rotor leaves a bivector alone, and a
    // rotor of norm k scales it by k^2.
    const Bivector<double> b(3.0);
    check_bivector(ga::sandwich(b, r), 3.0, "sandwich(B, unit r) == B");
    check_bivector(ga::rotate(b, r), 3.0, "rotate(B, r) == B");
    check_bivector(ga::sandwich(b, r * Scalar<double>(2.0)), 12.0, "sandwich(B, 2r) == 4B");

    // A multivector rotates grade by grade.
    const Multivector<double> m = make_mv(1.5, v.x, v.y, -2.0);
    const Multivector<double> rm = ga::sandwich(m, r);
    check_scalar(rm.scalar, 1.5, "sandwich(m, r) keeps the scalar");
    check_vector_close(rm.vector, ga::sandwich(v, r), "sandwich(m, r) rotates the vector");
    check_bivector(rm.bivector, -2.0, "sandwich(m, r) keeps the bivector");
    check_mv_close(ga::rotate(m, r), rm, "rotate(m, r) == sandwich(m, r)");
}

void test_slerp()
{
    section("slerp");

    const Rotor<double> a = ga::rotor_from_angle(0.2);
    const Rotor<double> b = ga::rotor_from_angle(1.0);

    check_rotor_close(ga::slerp(a, b, 0.0), a, "slerp at t = 0");
    check_rotor_close(ga::slerp(a, b, 1.0), b, "slerp at t = 1");
    check_close(ga::rotor_angle(ga::slerp(a, b, 0.25)), 0.4, "slerp interpolates the angle");
    check_close(ga::rotor_angle(ga::slerp(a, b, 0.5)), 0.6, "slerp midpoint");

    // -b is the same rotation as b; slerp must not take the long way round.
    check_close(ga::rotor_angle(ga::slerp(a, b * Scalar<double>(-1.0), 0.25)), 0.4,
                "slerp to -b takes the short arc");

    // Across the +-pi seam: from -3.0 to +3.0 the short way is through pi,
    // a turn of 2 pi - 6, not the 6 radians the other way.
    const Rotor<double> mid = ga::slerp(ga::rotor_from_angle(-3.0), ga::rotor_from_angle(3.0), 0.5);
    check_close(std::cos(ga::rotor_angle(mid)), -1.0, "slerp across the seam goes through pi", 1e-12);
    check_close(ga::norm(mid).value, 1.0, "slerp result is unit");
}

} // namespace

void run_rotations()
{
    test_exp_log();
    test_rotor_construction();
    test_sandwich();
    test_slerp();
}

} // namespace cctest
