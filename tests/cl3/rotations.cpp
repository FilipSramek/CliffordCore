// Cl(3,0) tests -- rotations: exp and log, rotation by the sandwich product on
// every grade, rotor construction and slerp.
//
// One category of the Cl3 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

void test_exp_log()
{
    section("exp/log");

    // exp of the zero bivector is the identity rotor.
    const Rotor<double> identity = ga::exp(Bivector<double>(0.0, 0.0, 0.0));
    check_scalar(identity.scalar, 1.0, "exp(0) scalar");
    check_bivector(identity.bivector, 0.0, 0.0, 0.0, "exp(0) bivector");

    // A quarter turn in the xy plane: cos(pi/2) + sin(pi/2) * xy.
    const Rotor<double> quarter = ga::exp(Bivector<double>(kPi / 2.0, 0.0, 0.0));
    check_close(quarter.scalar.value, 0.0, "exp(pi/2 xy) scalar", 1e-15);
    check_bivector(quarter.bivector, 1.0, 0.0, 0.0, "exp(pi/2 xy) bivector");

    // exp must always produce a unit rotor: scalar^2 + |bivector|^2 == 1.
    const Bivector<double> b(0.3, 0.4, 0.5);
    const Rotor<double> r = ga::exp(b);
    check_close(r.scalar.value * r.scalar.value + ga::squared_norm(r.bivector).value,
                1.0, "exp produces a unit rotor", 1e-15);

    // log of the identity rotor is the zero bivector.
    const Bivector<double> logIdentity =
        ga::log(Rotor<double>(Scalar<double>(1.0), Bivector<double>(0.0, 0.0, 0.0)));
    check_bivector(logIdentity, 0.0, 0.0, 0.0, "log(identity) is zero");

    // log and exp are inverses. Trig costs precision, so allow a looser bound.
    const Bivector<double> there = ga::log(ga::exp(b));
    check_bivector(there, b.xy, b.xz, b.yz, "log(exp(b)) == b");

    const Rotor<double> back = ga::exp(ga::log(r));
    check_close(back.scalar.value, r.scalar.value, "exp(log(r)) scalar", 1e-14);
    check_bivector(back.bivector, r.bivector.xy, r.bivector.xz, r.bivector.yz,
                   "exp(log(r)) bivector");

    // A rotor whose scalar has drifted just outside [-1, 1] must not yield NaN:
    // std::acos would be out of domain without the clamp in log().
    const Bivector<double> drifted =
        ga::log(Rotor<double>(Scalar<double>(1.0 + 2e-16), Bivector<double>(1.0, 0.0, 0.0)));
    check(!std::isnan(drifted.xy) && !std::isnan(drifted.xz) && !std::isnan(drifted.yz),
          "log clamps scalar drift above 1 instead of returning NaN");

    const Bivector<double> driftedLow =
        ga::log(Rotor<double>(Scalar<double>(-1.0 - 2e-16), Bivector<double>(1.0, 0.0, 0.0)));
    check(!std::isnan(driftedLow.xy), "log clamps scalar drift below -1");
}

void test_rotation()
{
    section("rotation");

    // exp(theta * B) rotates by 2*theta in the plane B, so pi/4 gives a quarter
    // turn. e12 spans the xy plane, so it turns about z and must leave e3 alone.
    const Rotor<double> quarterXY = ga::exp(Bivector<double>(kPi / 4.0, 0.0, 0.0));
    check_vector(ga::rotate(Vector<double>(1, 0, 0), quarterXY), 0.0, -1.0, 0.0,
                 "e1 rotates to -e2 in the xy plane", 1e-15);
    check_vector(ga::rotate(Vector<double>(0, 0, 1), quarterXY), 0.0, 0.0, 1.0,
                 "e3 is the axis of an xy rotation and is fixed", 1e-15);

    // e23 spans the yz plane, so it turns about x and must leave e1 alone.
    const Rotor<double> quarterYZ = ga::exp(Bivector<double>(0.0, 0.0, kPi / 4.0));
    check_vector(ga::rotate(Vector<double>(1, 0, 0), quarterYZ), 1.0, 0.0, 0.0,
                 "e1 is the axis of a yz rotation and is fixed", 1e-15);

    // Rotation preserves length.
    const Vector<double> v(1.0, 2.0, 3.0);
    check_scalar(ga::norm(ga::rotate(v, quarterXY)), std::sqrt(14.0),
                 "rotation preserves length");

    // The identity rotor leaves everything alone.
    const Rotor<double> identity(Scalar<double>(1.0), Bivector<double>(0.0, 0.0, 0.0));
    check_vector(ga::rotate(v, identity), v.x, v.y, v.z, "identity rotor fixes v");

    // The closed form must agree with the general product R v reverse(R).
    const Rotor<double> r =
        ga::normalize(Rotor<double>(Scalar<double>(0.3), Bivector<double>(0.5, -0.7, 0.2)));
    const Multivector<double> viaGeneral =
        rotor_to_mv(r) * make_mv(0, v.x, v.y, v.z, 0, 0, 0, 0) * rotor_to_mv(ga::reverse(r));
    const Vector<double> viaClosedForm = ga::sandwich(v, r);
    check_vector(viaClosedForm, viaGeneral.vector.x, viaGeneral.vector.y, viaGeneral.vector.z,
                 "sandwich matches R v reverse(R)", 1e-14);
    // The general product must come out a pure vector.
    check_close(viaGeneral.scalar.value, 0.0, "sandwich leaves no scalar part", 1e-14);
    check_close(viaGeneral.trivector.e123, 0.0, "sandwich leaves no trivector part", 1e-14);

    // Rotating twice by R equals rotating once by R*R.
    const Vector<double> twice = ga::rotate(ga::rotate(v, quarterXY), quarterXY);
    const Vector<double> once = ga::rotate(v, quarterXY * quarterXY);
    check_vector(twice, once.x, once.y, once.z, "rotating twice == rotating by the composed rotor", 1e-14);

    // A rotor and its inverse undo each other.
    const Vector<double> roundTrip =
        ga::rotate(ga::rotate(v, quarterXY), ga::inverse(quarterXY));
    check_vector(roundTrip, v.x, v.y, v.z, "rotate then rotate by the inverse returns v", 1e-14);
}

void test_rotor_construction()
{
    section("rotor construction");

    check_scalar(ga::identity_rotor<double>().scalar, 1.0, "identity rotor scalar");
    check_bivector(ga::identity_rotor<double>().bivector, 0.0, 0.0, 0.0,
                   "identity rotor bivector");

    // Right-hand rule: +90 degrees about +z takes e1 to e2.
    const Rotor<double> aboutZ =
        ga::rotor_from_axis_angle(Vector<double>(0, 0, 1), kPi / 2.0);
    check_vector(ga::rotate(Vector<double>(1, 0, 0), aboutZ), 0.0, 1.0, 0.0,
                 "+90 about z takes e1 to e2", 1e-15);
    check_vector(ga::rotate(Vector<double>(0, 1, 0), aboutZ), -1.0, 0.0, 0.0,
                 "+90 about z takes e2 to -e1", 1e-15);
    check_vector(ga::rotate(Vector<double>(0, 0, 1), aboutZ), 0.0, 0.0, 1.0,
                 "the axis is fixed", 1e-15);

    // +90 about +x takes e2 to e3.
    const Rotor<double> aboutX =
        ga::rotor_from_axis_angle(Vector<double>(1, 0, 0), kPi / 2.0);
    check_vector(ga::rotate(Vector<double>(0, 1, 0), aboutX), 0.0, 0.0, 1.0,
                 "+90 about x takes e2 to e3", 1e-15);

    // The axis need not be unit length.
    const Rotor<double> unnormalised =
        ga::rotor_from_axis_angle(Vector<double>(0, 0, 5), kPi / 2.0);
    check_vector(ga::rotate(Vector<double>(1, 0, 0), unnormalised), 0.0, 1.0, 0.0,
                 "axis is normalized internally", 1e-15);

    // A zero angle is the identity; a zero axis degrades to the identity too.
    check_vector(ga::rotate(Vector<double>(1, 2, 3),
                     ga::rotor_from_axis_angle(Vector<double>(0, 0, 1), 0.0)),
                 1.0, 2.0, 3.0, "zero angle is the identity");
    check_scalar(ga::rotor_from_axis_angle(Vector<double>(0, 0, 0), 1.0).scalar, 1.0,
                 "zero axis degrades to the identity");

    // rotor_between must actually carry one direction onto the other.
    const Vector<double> from(1, 2, 3);
    const Vector<double> to(-2, 1, 0.5);
    const Rotor<double> between = ga::rotor_between(from, to);
    const Vector<double> carried = ga::rotate(ga::normalize(from), between);
    const Vector<double> target = ga::normalize(to);
    check_vector(carried, target.x, target.y, target.z, "rotor_between carries from onto to", 1e-14);
    check_scalar(ga::norm(between), 1.0, "rotor_between returns a unit rotor");

    // Identical directions give no rotation.
    check_vector(ga::rotate(Vector<double>(1, 0, 0),
                     ga::rotor_between(Vector<double>(1, 0, 0), Vector<double>(2, 0, 0))),
                 1.0, 0.0, 0.0, "rotor_between of parallel vectors is the identity", 1e-14);

    // Opposite directions are the degenerate case: the plane is ambiguous, but
    // the result must still be a genuine half turn rather than NaN.
    const Rotor<double> flip =
        ga::rotor_between(Vector<double>(1, 0, 0), Vector<double>(-1, 0, 0));
    check_vector(ga::rotate(Vector<double>(1, 0, 0), flip), -1.0, 0.0, 0.0,
                 "antiparallel gives a half turn", 1e-14);
    check_scalar(ga::norm(flip), 1.0, "the half turn is a unit rotor");

    // slerp interpolates: the endpoints are exact and the midpoint is halfway.
    const Rotor<double> identity = ga::identity_rotor<double>();
    check_vector(ga::rotate(Vector<double>(1, 0, 0),
                     ga::slerp(identity, aboutZ, 0.0)),
                 1.0, 0.0, 0.0, "slerp at t=0 is the start", 1e-14);
    check_vector(ga::rotate(Vector<double>(1, 0, 0),
                     ga::slerp(identity, aboutZ, 1.0)),
                 0.0, 1.0, 0.0, "slerp at t=1 is the end", 1e-14);
    const double half = std::sqrt(0.5);
    check_vector(ga::rotate(Vector<double>(1, 0, 0),
                     ga::slerp(identity, aboutZ, 0.5)),
                 half, half, 0.0, "slerp at t=0.5 is a 45 degree turn", 1e-14);
}

void test_sandwich_all_grades()
{
    section("sandwich (all grades)");

    const Rotor<double> r = ga::rotor_from_axis_angle(Vector<double>(0, 0, 1), kPi / 2);

    // Rotating a plane. The e13 plane, turned a quarter turn about z, is e23.
    check_bivector(ga::rotate(Bivector<double>(0, 1, 0), r), 0.0, 0.0, 1.0,
                   "the e13 plane rotates to e23", 1e-15);

    // The rotation plane itself is fixed, exactly as the axis vector is.
    check_bivector(ga::rotate(Bivector<double>(1, 0, 0), r), 1.0, 0.0, 0.0,
                   "the e12 plane is the rotation plane and is fixed", 1e-15);

    // The pseudoscalar commutes with everything, so it cannot be rotated.
    check_trivector(ga::rotate(Trivector<double>(5), r), 5.0,
                    "the pseudoscalar is central and survives untouched", 1e-14);

    // Rotating a multivector must agree grade by grade with rotating the parts.
    const Multivector<double> m(Scalar<double>(1), Vector<double>(1, 0, 0),
                                 Bivector<double>(0, 1, 0), Trivector<double>(5));
    const Multivector<double> rotated = ga::rotate(m, r);
    check_scalar(rotated.scalar, 1.0, "the scalar part is untouched");
    const Vector<double> rv = ga::rotate(Vector<double>(1, 0, 0), r);
    check_vector(rotated.vector, rv.x, rv.y, rv.z, "grade 1 matches rotating the vector alone", 1e-15);
    const Bivector<double> rb = ga::rotate(Bivector<double>(0, 1, 0), r);
    check_bivector(rotated.bivector, rb.xy, rb.xz, rb.yz,
                   "grade 2 matches rotating the bivector alone", 1e-15);

    // The four operator* pairs that were missing, without which none of the
    // above could be written.
    const Bivector<double> b(1, 2, 3);
    const Trivector<double> t(4);
    check_close(mv_difference(r * b, ga::to_multivector(r) * b), 0.0, "Rotor * Bivector");
    check_close(mv_difference(b * r, b * ga::to_multivector(r)), 0.0, "Bivector * Rotor");
    check_close(mv_difference(r * t, ga::to_multivector(r) * t), 0.0, "Rotor * Trivector");
    check_close(mv_difference(t * r, t * ga::to_multivector(r)), 0.0, "Trivector * Rotor");
}

} // namespace

void run_rotations()
{
    test_exp_log();
    test_rotation();
    test_rotor_construction();
    test_sandwich_all_grades();
}

} // namespace cctest
