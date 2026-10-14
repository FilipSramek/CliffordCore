// Cl(2,0) tests -- geometry: reflect, project and reject, and the planar
// extras: perp, signed_area, angle and angle_between.
//
// One category of the Cl2 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

void test_geometry()
{
    section("geometry");

    const Vector<double> e1(1, 0), e2(0, 1);
    const Vector<double> v(3.0, 1.0);
    const Vector<double> diagonal(1.0, 1.0);

    // reflect(v, n) mirrors in the line perpendicular to n. Values from the
    // independent model.
    check_vector(ga::reflect(e1, e1), -1.0, 0.0, "reflect(e1, n = e1)");
    check_vector(ga::reflect(e2, e1), 0.0, 1.0, "reflect(e2, n = e1)");
    check_vector(ga::reflect(v, e1), -3.0, 1.0, "reflect((3,1), n = e1)");
    check_vector(ga::reflect(v, e2), 3.0, -1.0, "reflect((3,1), n = e2)");
    check_vector(ga::reflect(e1, diagonal), 0.0, -1.0, "reflect(e1, n = (1,1))");
    check_vector(ga::reflect(v, diagonal), -1.0, -3.0, "reflect((3,1), n = (1,1))");
    check_vector_close(ga::reflect(v, diagonal * Scalar<double>(5.0)), ga::reflect(v, diagonal),
                       "the normal need not be unit");
    check_vector_close(ga::reflect(ga::reflect(v, diagonal), diagonal), v, "reflecting twice is the identity");

    // Reflecting in the line ALONG u is -reflect(v, u).
    check_vector(-ga::reflect(v, e1), 3.0, -1.0, "-reflect(v, e1) mirrors in the x axis");

    // Two reflections make a rotation, by twice the angle from the first
    // normal to the second.
    const Vector<double> n2(std::cos(0.3), std::sin(0.3));
    check_vector_close(ga::reflect(ga::reflect(v, e1), n2), ga::rotate(v, ga::rotor_from_angle(0.6)),
                       "two reflections == a rotation by twice the angle");

    // Projection and rejection split a vector exactly.
    const Vector<double> u(2.0, 0.0);
    check_vector(ga::project(Vector<double>(3, 4), u), 3.0, 0.0, "project onto the x axis");
    check_vector(ga::reject(Vector<double>(3, 4), u), 0.0, 4.0, "reject from the x axis");
    const Vector<double> t(1.0, -2.0);
    check_vector_close(ga::project(v, t) + ga::reject(v, t), v, "project + reject == v");
    check_close((ga::reject(v, t) | t).value, 0.0, "reject is perpendicular to u");
    check_close((ga::project(v, t) ^ t).xy, 0.0, "project is parallel to u");

    // perp is the counter-clockwise quarter turn, and it is exactly dual().
    check_vector(ga::perp(Vector<double>(3, 4)), -4.0, 3.0, "perp((3,4)) == (-4,3)");
    check_vector_close(ga::perp(v), ga::dual(v), "perp == dual");
    check_close((ga::perp(v) | v).value, 0.0, "perp(v) is perpendicular to v");
    check_vector_close(ga::perp(ga::perp(v)), -v, "perp twice negates");
    check_vector_close(ga::perp(v), ga::rotate(v, ga::rotor_from_angle(kPi / 2)), "perp is +90 degrees");

    // signed_area is positive counter-clockwise, and is the wedge component.
    check_scalar(ga::signed_area(e1, e2), 1.0, "signed_area(e1, e2) == +1");
    check_scalar(ga::signed_area(e2, e1), -1.0, "signed_area(e2, e1) == -1");
    check_scalar(ga::signed_area(Vector<double>(1, 2), Vector<double>(3, 4)), -2.0, "signed_area((1,2),(3,4))");
    check_scalar(ga::signed_area(v, v), 0.0, "signed_area(v, v) == 0");
    check_scalar(ga::signed_area(v, ga::perp(v)), ga::squared_norm(v).value, "signed_area(v, perp v) == |v|^2");

    // angle_between is signed; angle is the direction from +x.
    check_close(ga::angle_between(e1, e2), kPi / 2, "angle_between(e1, e2) == +pi/2");
    check_close(ga::angle_between(e2, e1), -kPi / 2, "angle_between(e2, e1) == -pi/2");
    check_close(ga::angle_between(e1, -e1), kPi, "angle_between(e1, -e1) == pi");
    check_close(ga::angle_between(e1 * Scalar<double>(4.0), diagonal), kPi / 4, "angle_between ignores length");
    check_close(ga::angle(e1), 0.0, "angle(e1) == 0");
    check_close(ga::angle(e2), kPi / 2, "angle(e2) == pi/2");
    check_close(ga::angle(-e1), kPi, "angle(-e1) == pi");
    check_close(ga::angle(-e2), -kPi / 2, "angle(-e2) == -pi/2");

    // The principal range is (-pi, pi], and signed zero is what threatens it:
    // each input below makes plain atan2 return -pi. -e1 is (-1, -0.0), and
    // (-1, 0) ^ (1, 0) is -0.0.
    check(ga::angle(-e1) > 0.0, "angle(-e1) is +pi, not -pi");
    check(ga::angle_between(Vector<double>(-1, 0), e1) > 0.0, "a half turn between vectors is +pi");
    check_close(ga::angle_between(Vector<double>(-1, 0), e1), kPi, "angle_between((-1,0), e1) == pi");
}

} // namespace

void run_geometry()
{
    test_geometry();
}

} // namespace cctest
