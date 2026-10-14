// Cl(3,0,1) tests -- motions: translator, rotor and motor construction, exp
// and log, and slerp for all three versor types.
//
// One category of the PGA suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

// A unit motor satisfies m * reverse(m) == 1 exactly -- scalar one and nothing
// else, e0123 included. That last part is what a naive normalisation misses.
void check_motor_unit(const Motor<double>& m, const std::string& what, double tolerance = 1e-12)
{
    const Motor<double> product = m * ga::reverse(m);
    check_close(product.scalar.value, 1.0, what + " (scalar)", tolerance);
    check_close(product.quadvector.e0123, 0.0, what + " (e0123)", tolerance);
}

void test_translator_construction()
{
    section("translator");

    // The defining behaviour: it moves a point by exactly (dx, dy, dz).
    const Trivector<double> O = ga::origin<double>();
    check_trivector(ga::sandwich(O, ga::translator(2.0, 0.0, 0.0)), ga::point(2.0, 0.0, 0.0).e012,
                    ga::point(2.0, 0.0, 0.0).e013, ga::point(2.0, 0.0, 0.0).e023, 1.0, "translator moves the origin along x");
    check_trivector(ga::sandwich(O, ga::translator(0.0, 2.0, 0.0)), 0, 2, 0, 1, "translator moves the origin along y");
    check_trivector(ga::sandwich(O, ga::translator(0.0, 0.0, 2.0)), -2, 0, 0, 1, "translator moves the origin along z");
    const Trivector<double> moved = ga::sandwich(ga::point(1.0, 1.0, 1.0), ga::translator(1.0, 2.0, 3.0));
    check_trivector(moved, ga::point(2.0, 3.0, 4.0).e012, ga::point(2.0, 3.0, 4.0).e013,
                    ga::point(2.0, 3.0, 4.0).e023, 1.0, "translator moves an off-origin point");

    // Planes and lines move too.
    check_vector(ga::sandwich(ga::plane(1.0, 0.0, 0.0, 0.0), ga::translator(1.0, 0.0, 0.0)), -1, 1, 0, 0,
                 "translating x = 0 by +1 gives x = 1");
    const Bivector<double> zAxis = ga::line_through_origin(0.0, 0.0, 1.0);
    check_bivector(ga::translate(zAxis, ga::translator(1.0, 0.0, 0.0)), 0, -1, 0, 1, 0, 0,
                   "translating the z axis by +x gives the z-parallel line through (1,0,0)");
    check_close(ga::distance(ga::point(0.0, 0.0, 5.0), ga::translate(zAxis, ga::translator(1.0, 0.0, 0.0))), 1.0,
                "... which is 1 away from the z axis");

    // Composition, inversion, and the named constructors.
    const Translator<double> a = ga::translator(1.0, 2.0, 3.0);
    const Translator<double> b = ga::translator(-1.0, 4.0, 0.5);
    check_translator(a * b, ga::translator(0.0, 6.0, 3.5).scalar.value, ga::translator(0.0, 6.0, 3.5).e01,
                     ga::translator(0.0, 6.0, 3.5).e02, ga::translator(0.0, 6.0, 3.5).e03, "translators add");
    check_trivector(ga::sandwich(ga::sandwich(O, a), ga::inverse(a)), 0, 0, 0, 1, "a translation undoes itself");
    check_translator(ga::identity_translator<double>(), 1, 0, 0, 0, "identity_translator");
    check_trivector(ga::sandwich(ga::point(3.0, -1.0, 2.0), ga::identity_translator<double>()), -2, -1, -3, 1,
                    "the identity translator moves nothing");
    const Translator<double> between = ga::translator_between(ga::point(1.0, 2.0, 3.0), ga::point(4.0, 6.0, 3.0));
    check_close(ga::distance(ga::sandwich(ga::point(1.0, 2.0, 3.0), between), ga::point(4.0, 6.0, 3.0)), 0.0,
                "translator_between lands exactly");
}

void test_rotor_construction()
{
    section("rotor construction");

    const double q = kPi / 2;

    // The right-hand rule, on all three axes, for points.
    check_close(ga::distance(ga::sandwich(ga::point(1.0, 0.0, 0.0), ga::rotor_from_axis_angle(0.0, 0.0, 1.0, q)),
                             ga::point(0.0, 1.0, 0.0)), 0.0, "+90 about +z takes (1,0,0) to (0,1,0)");
    check_close(ga::distance(ga::sandwich(ga::point(0.0, 1.0, 0.0), ga::rotor_from_axis_angle(1.0, 0.0, 0.0, q)),
                             ga::point(0.0, 0.0, 1.0)), 0.0, "+90 about +x takes (0,1,0) to (0,0,1)");
    check_close(ga::distance(ga::sandwich(ga::point(0.0, 0.0, 1.0), ga::rotor_from_axis_angle(0.0, 1.0, 0.0, q)),
                             ga::point(1.0, 0.0, 0.0)), 0.0, "+90 about +y takes (0,0,1) to (1,0,0)");

    // ... and for planes, which is the same rotation seen the other way up.
    const Vector<double> turned = ga::rotate(ga::plane(1.0, 0.0, 0.0, 0.0), ga::rotor_from_axis_angle(0.0, 0.0, 1.0, q));
    check_vector(turned, 0, 0, 1, 0, "+90 about +z takes the plane x=0 to the plane y=0");

    // The axis is fixed, and the rotor is unit.
    const Rotor<double> r = ga::rotor_from_axis_angle(1.0, 1.0, 1.0, 0.7);
    check_scalar(ga::norm(r), 1.0, "rotor_from_axis_angle returns a unit rotor");
    const Bivector<double> axis = ga::normalize(ga::line_through_origin(1.0, 1.0, 1.0));
    check_bivector(ga::rotate(axis, r), axis.e01, axis.e02, axis.e03, axis.e12, axis.e13, axis.e23,
                   "the axis line is fixed by the rotation");
    check_trivector(ga::sandwich(ga::origin<double>(), r), 0, 0, 0, 1, "the origin is fixed by a rotor");

    // An unnormalised axis, a zero axis, and composition.
    check_rotor(ga::rotor_from_axis_angle(0.0, 0.0, 5.0, q), ga::rotor_from_axis_angle(0.0, 0.0, 1.0, q).scalar.value,
                ga::rotor_from_axis_angle(0.0, 0.0, 1.0, q).e12, ga::rotor_from_axis_angle(0.0, 0.0, 1.0, q).e13,
                ga::rotor_from_axis_angle(0.0, 0.0, 1.0, q).e23, "the axis is normalized internally");
    check_rotor(ga::rotor_from_axis_angle(0.0, 0.0, 0.0, q), 1, 0, 0, 0, "a zero axis gives the identity");
    check_rotor(ga::identity_rotor<double>(), 1, 0, 0, 0, "identity_rotor");
    const Rotor<double> half = ga::rotor_from_axis_angle(0.0, 0.0, 1.0, q / 2);
    check_close(ga::distance(ga::sandwich(ga::point(1.0, 0.0, 0.0), half * half), ga::point(0.0, 1.0, 0.0)), 0.0,
                "two 45 degree turns compose into 90");
    check_close(ga::distance(ga::sandwich(ga::sandwich(ga::point(1.0, 2.0, 3.0), r), ga::reverse(r)), ga::point(1.0, 2.0, 3.0)),
                0.0, "reversing a unit rotor undoes it");

    // A rotor agrees with the motor built from the same line.
    const Motor<double> asMotor = ga::motor_from_line_angle(ga::line_through_origin(0.0, 0.0, 1.0), q);
    check_mv_close(ga::to_multivector(asMotor), ga::to_multivector(Motor<double>(ga::rotor_from_axis_angle(0.0, 0.0, 1.0, q))),
                   "motor_from_line_angle about an origin line is the rotor");
}

void test_motor_construction()
{
    section("motor construction");

    const double q = kPi / 2;

    // A rotation about a line that misses the origin: the whole point of PGA.
    const Bivector<double> offset = ga::line_through_points(ga::point(1.0, 0.0, 0.0), ga::point(1.0, 0.0, 1.0));
    const Motor<double> turn = ga::motor_from_line_angle(offset, kPi);
    check_close(ga::distance(ga::sandwich(ga::point(2.0, 0.0, 0.0), turn), ga::origin<double>()), 0.0,
                "180 degrees about the z-parallel line through (1,0,0) maps (2,0,0) to the origin");
    check_close(ga::distance(ga::sandwich(ga::point(1.0, 0.0, 5.0), turn), ga::point(1.0, 0.0, 5.0)), 0.0,
                "... and fixes the points on that line");
    check_scalar(ga::norm(turn), 1.0, "motor_from_line_angle returns a unit motor");
    check_motor_unit(turn, "motor_from_line_angle: m * reverse(m) == 1");

    // A screw: rotate about the z axis and rise along it at the same time.
    const Motor<double> s = ga::screw(ga::line_through_origin(0.0, 0.0, 1.0), q, 2.0);
    const Trivector<double> screwed = ga::sandwich(ga::point(1.0, 0.0, 0.0), s);
    check_close(ga::distance(screwed, ga::point(0.0, 1.0, 2.0)), 0.0, "screw turns by 90 about z and rises by 2");
    check_motor_unit(s, "screw: m * reverse(m) == 1");
    check_close(ga::distance(ga::sandwich(ga::origin<double>(), ga::screw(ga::line_through_origin(0.0, 0.0, 1.0), 0.0, 3.0)),
                             ga::point(0.0, 0.0, 3.0)), 0.0, "a screw with no rotation is a translation");
    check_mv_close(ga::to_multivector(ga::screw(ga::line_through_origin(0.0, 0.0, 1.0), q, 0.0)),
                   ga::to_multivector(ga::motor_from_line_angle(ga::line_through_origin(0.0, 0.0, 1.0), q)),
                   "a screw with no rise is a rotation");
    check_mv_close(ga::to_multivector(ga::screw(Bivector<double>(1, 2, 3, 0, 0, 0), q, 1.0)),
                   ga::to_multivector(ga::identity_motor<double>()), "an ideal axis gives the identity motor");

    // motor_between takes a plane onto a plane.
    const Vector<double> from = ga::plane(1.0, 0.0, 0.0, 0.0);
    const Vector<double> to = ga::plane(0.0, 1.0, 0.0, 0.0);
    const Motor<double> between = ga::motor_between(from, to);
    const Vector<double> landed = ga::sandwich(from, between);
    check_vector(landed, to.e0, to.e1, to.e2, to.e3, "motor_between lands the plane exactly");
    check_motor_unit(between, "motor_between returns a unit motor");
    const Motor<double> opposite = ga::motor_between(from, -from);
    check_vector(ga::sandwich(from, opposite), -from.e0, -from.e1, -from.e2, -from.e3,
                 "motor_between handles opposite planes with a half turn");
    // Opposite normals at different offsets: the half turn has to be about a
    // line halfway between the planes. A line in `from` -- what this once used --
    // lands on -from, x = 0 facing -x, instead of on x = 2.
    const Vector<double> facing_back = ga::plane(-1.0, 0.0, 0.0, 2.0);   // x = 2, facing -x
    check_vector(ga::sandwich(from, ga::motor_between(from, facing_back)), 2, -1, 0, 0,
                 "motor_between handles opposite planes at different offsets");
    check_mv_close(ga::to_multivector(ga::motor_between(from, ga::plane(0.0, 0.0, 0.0, 0.0))),
                   ga::to_multivector(ga::identity_motor<double>()), "a degenerate target gives the identity");

    // Rotation then translation, and the identity.
    const Motor<double> rt = ga::motor_from_rotor_translator(ga::rotor_from_axis_angle(0.0, 0.0, 1.0, q),
                                                             ga::translator(0.0, 0.0, 5.0));
    check_close(ga::distance(ga::sandwich(ga::point(1.0, 0.0, 0.0), rt), ga::point(0.0, 1.0, 5.0)), 0.0,
                "motor_from_rotor_translator rotates then translates");
    check_close(ga::distance(ga::sandwich(ga::point(1.0, 2.0, 3.0), ga::identity_motor<double>()), ga::point(1.0, 2.0, 3.0)),
                0.0, "the identity motor moves nothing");

    // Motors compose and invert.
    const Motor<double> a = ga::screw(ga::line_through_points(ga::point(1.0, 0.0, 0.0), ga::point(1.0, 1.0, 0.0)), 0.6, 1.5);
    const Motor<double> b = ga::motor_from_line_angle(ga::line_through_origin(1.0, 1.0, 1.0), -0.9);
    const Trivector<double> p = ga::point(2.0, -1.0, 0.5);
    check_close(ga::distance(ga::sandwich(p, a * b), ga::sandwich(ga::sandwich(p, b), a)), 0.0,
                "(a*b) applied once == b applied then a");
    check_close(ga::distance(ga::sandwich(ga::sandwich(p, a), ga::inverse(a)), p), 0.0, "a motion undoes itself");
    check_close(ga::distance(ga::sandwich(ga::sandwich(p, a), ga::reverse(a)), p), 0.0, "reverse undoes a unit motor");

    // A rigid motion preserves distances and angles.
    const Trivector<double> p2 = ga::point(-1.0, 3.0, 2.0);
    check_close(ga::distance(ga::sandwich(p, a), ga::sandwich(p2, a)), ga::distance(p, p2), "a motor preserves distances");
    check_close(ga::angle(ga::sandwich(from, a), ga::sandwich(to, a)), ga::angle(from, to), "a motor preserves angles");
}

void test_exp_log()
{
    section("exp/log");

    // exp of zero, of a pure Euclidean bivector, and of a pure ideal one.
    check_mv_close(ga::to_multivector(ga::exp(Bivector<double>())), ga::to_multivector(ga::identity_motor<double>()),
                   "exp(0) is the identity motor");

    const Bivector<double> euclid(0, 0, 0, 0.7, -0.2, 0.4);
    const double theta = std::sqrt(0.49 + 0.04 + 0.16);
    const Motor<double> fromEuclid = ga::exp(euclid);
    check_scalar(fromEuclid.scalar, std::cos(theta), "exp of a Euclidean bivector: scalar is cos|B|");
    check_bivector(fromEuclid.bivector, 0, 0, 0, std::sin(theta) * 0.7 / theta, std::sin(theta) * -0.2 / theta,
                   std::sin(theta) * 0.4 / theta, "... bivector is sin|B| B/|B|, as in Cl3");
    check_quadvector(fromEuclid.quadvector, 0.0, "... and there is no e0123 part");

    const Bivector<double> ideal(1.5, -2.0, 0.5, 0, 0, 0);
    const Motor<double> fromIdeal = ga::exp(ideal);
    check_scalar(fromIdeal.scalar, 1.0, "exp of an ideal bivector: scalar is exactly 1");
    check_bivector(fromIdeal.bivector, 1.5, -2.0, 0.5, 0, 0, 0, "... and the series terminates at 1 + B");
    check_quadvector(fromIdeal.quadvector, 0.0, "... with no e0123 part");

    // A general screw: exp is unit, and log inverts it.
    const Bivector<double> screwB(0.3, -0.7, 1.1, 0.5, 0.2, -0.9);
    const Motor<double> m = ga::exp(screwB);
    check_motor_unit(m, "exp(B) is always a unit motor");
    check_bivector(ga::log(m), screwB.e01, screwB.e02, screwB.e03, screwB.e12, screwB.e13, screwB.e23,
                   "log(exp(B)) == B for a general screw", 1e-12);
    check_mv_close(ga::to_multivector(ga::exp(ga::log(m))), ga::to_multivector(m), "exp(log(m)) == m");
    check_bivector(ga::log(ga::exp(euclid)), 0, 0, 0, 0.7, -0.2, 0.4, "log(exp(B)) == B for a Euclidean bivector");
    check_bivector(ga::log(ga::exp(ideal)), 1.5, -2.0, 0.5, 0, 0, 0, "log(exp(B)) == B for an ideal bivector");

    // exp and the constructors agree.
    const double q = kPi / 2;
    const Bivector<double> zAxis = ga::line_through_origin(0.0, 0.0, 1.0);
    check_mv_close(ga::to_multivector(ga::exp(zAxis * Scalar<double>(-q / 2))),
                   ga::to_multivector(ga::motor_from_line_angle(zAxis, q)), "exp(-(angle/2) L) == motor_from_line_angle");
    check_bivector(ga::log(ga::translator(2.0, 4.0, 6.0)), -1, -2, -3, 0, 0, 0, "log(translator) is its ideal part");
    // exp(-(angle/2) L) means the log of a quarter turn about z is -(pi/4) e12.
    check_bivector(ga::log(ga::rotor_from_axis_angle(0.0, 0.0, 1.0, q)), 0, 0, 0, -q / 2, 0, 0,
                   "log(rotor) recovers half the angle on the rotation plane");

    // log recovers the axis and angle of a motion.
    const Motor<double> known = ga::motor_from_line_angle(ga::normalize(ga::line_through_origin(1.0, 2.0, 2.0)), 1.2);
    const Bivector<double> recovered = ga::log(known);
    check_close(2 * ga::norm(recovered).value, 1.2, "log recovers the angle, as twice the bivector norm");
    const std::array<double, 3> axisDir = ga::direction(ga::normalize(recovered));
    check_close(std::abs(axisDir[0]), 1.0 / 3.0, "log recovers the axis direction x");
    check_close(std::abs(axisDir[1]), 2.0 / 3.0, "... y");
    check_close(std::abs(axisDir[2]), 2.0 / 3.0, "... z");

    // Zero motors and the degenerate cases do not produce NaN.
    check_bivector(ga::log(Motor<double>()), 0, 0, 0, 0, 0, 0, "log of a zero motor returns zero, not NaN");
    check_bivector(ga::log(ga::identity_motor<double>()), 0, 0, 0, 0, 0, 0, "log(identity) is zero");
}

void test_slerp()
{
    section("slerp");

    const double q = kPi / 2;
    const Motor<double> a = ga::identity_motor<double>();
    const Motor<double> b = ga::motor_from_line_angle(ga::line_through_origin(0.0, 0.0, 1.0), q);

    check_mv_close(ga::to_multivector(ga::slerp(a, b, 0.0)), ga::to_multivector(a), "slerp at t = 0 is the start");
    check_mv_close(ga::to_multivector(ga::slerp(a, b, 1.0)), ga::to_multivector(b), "slerp at t = 1 is the end", 1e-12);
    check_close(ga::distance(ga::sandwich(ga::point(1.0, 0.0, 0.0), ga::slerp(a, b, 0.5)),
                             ga::sandwich(ga::point(1.0, 0.0, 0.0), ga::motor_from_line_angle(ga::line_through_origin(0.0, 0.0, 1.0), q / 2))),
                0.0, "slerp at t = 0.5 is the half rotation");
    check_motor_unit(ga::slerp(a, b, 0.37), "slerp stays on the unit motors");

    // A screw interpolates rotation and translation together.
    const Motor<double> s = ga::screw(ga::line_through_origin(0.0, 0.0, 1.0), q, 4.0);
    const Trivector<double> mid = ga::sandwich(ga::point(1.0, 0.0, 0.0), ga::slerp(a, s, 0.5));
    const std::array<double, 3> midPos = ga::position(mid);
    check_close(midPos[2], 2.0, "half a screw has risen half the distance");
    check_close(std::atan2(midPos[1], midPos[0]), q / 2, "... and turned half the angle");

    // The short way round: a negated motor is the same motion.
    check_close(ga::distance(ga::sandwich(ga::point(1.0, 0.0, 0.0), ga::slerp(a, -b, 0.5)),
                             ga::sandwich(ga::point(1.0, 0.0, 0.0), ga::slerp(a, b, 0.5))), 0.0,
                "slerp to a negated motor takes the same short path");

    // The rotor and translator overloads.
    const Rotor<double> r0 = ga::identity_rotor<double>();
    const Rotor<double> r1 = ga::rotor_from_axis_angle(0.0, 0.0, 1.0, q);
    check_rotor(ga::slerp(r0, r1, 0.0), 1, 0, 0, 0, "slerp(rotor) at t = 0");
    check_close(ga::distance(ga::sandwich(ga::point(1.0, 0.0, 0.0), ga::slerp(r0, r1, 1.0)), ga::point(0.0, 1.0, 0.0)), 0.0,
                "slerp(rotor) at t = 1");
    static_assert(std::is_same<decltype(ga::slerp(r0, r1, 0.5)), Rotor<double>>::value, "slerp(Rotor, Rotor) returns a Rotor");
    const Translator<double> t1 = ga::translator(4.0, 0.0, 0.0);
    check_close(ga::distance(ga::sandwich(ga::origin<double>(), ga::slerp(ga::identity_translator<double>(), t1, 0.25)),
                             ga::point(1.0, 0.0, 0.0)), 0.0, "slerp(translator) moves at constant speed");
    static_assert(std::is_same<decltype(ga::slerp(t1, t1, 0.5)), Translator<double>>::value,
                  "slerp(Translator, Translator) returns a Translator");
}

} // namespace

void run_motions()
{
    test_translator_construction();
    test_rotor_construction();
    test_motor_construction();
    test_exp_log();
    test_slerp();
}

} // namespace cctest
