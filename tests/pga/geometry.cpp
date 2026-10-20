// Cl(3,0,1) tests -- geometry: the primitives.hpp factories and accessors,
// then reflect, project, distance and angle.
//
// One category of the PGA suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

void test_primitives()
{
    section("primitives");

    // The conventions from docs/pga.md, spelled out once.
    check_trivector(ga::point(1.0, 2.0, 3.0), -3, 2, -1, 1, "point(1,2,3) == e123 - x e023 + y e013 - z e012");
    check_trivector(ga::origin<double>(), 0, 0, 0, 1, "origin() == e123");
    check_trivector(ga::ideal_point(1.0, 2.0, 3.0), -3, 2, -1, 0, "ideal_point has the same ideal parts and zero weight");
    check_vector(ga::plane(1.0, 2.0, 3.0, 4.0), 4, 1, 2, 3, "plane(a,b,c,d) == d e0 + a e1 + b e2 + c e3");
    check_vector(ga::ideal_plane<double>(), 1, 0, 0, 0, "ideal_plane() == e0");
    check_bivector(ga::line_through_origin(1.0, 0.0, 0.0), 0, 0, 0, 0, 0, 1, "the x axis is e23");
    check_bivector(ga::line_through_origin(0.0, 1.0, 0.0), 0, 0, 0, 0, -1, 0, "the y axis is -e13");
    check_bivector(ga::line_through_origin(0.0, 0.0, 1.0), 0, 0, 0, 1, 0, 0, "the z axis is e12");

    // A point built by the factory is the triple meet of its coordinate planes.
    check_trivector(ga::point(1.0, 2.0, 3.0),
                    (ga::plane(1.0, 0.0, 0.0, -1.0) ^ ga::plane(0.0, 1.0, 0.0, -2.0) ^ ga::plane(0.0, 0.0, 1.0, -3.0)).e012,
                    (ga::plane(1.0, 0.0, 0.0, -1.0) ^ ga::plane(0.0, 1.0, 0.0, -2.0) ^ ga::plane(0.0, 0.0, 1.0, -3.0)).e013,
                    (ga::plane(1.0, 0.0, 0.0, -1.0) ^ ga::plane(0.0, 1.0, 0.0, -2.0) ^ ga::plane(0.0, 0.0, 1.0, -3.0)).e023,
                    (ga::plane(1.0, 0.0, 0.0, -1.0) ^ ga::plane(0.0, 1.0, 0.0, -2.0) ^ ga::plane(0.0, 0.0, 1.0, -3.0)).e123,
                    "point(1,2,3) is the meet of x=1, y=2, z=3");

    // Round trips.
    const std::array<double, 3> pos = ga::position(ga::point(1.5, -2.5, 3.0));
    check_close(pos[0], 1.5, "position round trips x");
    check_close(pos[1], -2.5, "position round trips y");
    check_close(pos[2], 3.0, "position round trips z");
    const std::array<double, 3> scaled = ga::position(ga::point(1.5, -2.5, 3.0) * Scalar<double>(4.0));
    check_close(scaled[0], 1.5, "position divides out the weight");
    check_close(scaled[2], 3.0, "... on every coordinate");
    const std::array<double, 3> none = ga::position(ga::ideal_point(1.0, 2.0, 3.0));
    check(none[0] == 0.0 && none[1] == 0.0 && none[2] == 0.0, "an ideal point has no position");

    const std::array<double, 3> dir = ga::direction(ga::line_through_origin(1.0, 2.0, 3.0));
    check_close(dir[0], 1.0, "direction(line) round trips dx");
    check_close(dir[1], 2.0, "direction(line) round trips dy");
    check_close(dir[2], 3.0, "direction(line) round trips dz");
    const std::array<double, 3> nrm = ga::normal(ga::plane(1.0, 2.0, 3.0, 4.0));
    check_close(nrm[0], 1.0, "normal(plane) round trips a");
    check_close(nrm[2], 3.0, "normal(plane) round trips c");
    check_scalar(ga::offset(ga::plane(1.0, 2.0, 3.0, 4.0)), 4.0, "offset(plane) is d");
    const std::array<double, 3> mom = ga::moment(Bivector<double>(1, 2, 3, 4, 5, 6));
    check(mom[0] == 1.0 && mom[1] == 2.0 && mom[2] == 3.0, "moment(line) is the ideal half");

    // Lines two ways, and the plane through three points.
    check_bivector(ga::line_through_points(ga::origin<double>(), ga::point(0.0, 0.0, 1.0)), 0, 0, 0, 1, 0, 0,
                   "line_through_points(O, Z) is the z axis");
    check_bivector(ga::line_from_planes(ga::plane(1.0, 0.0, 0.0, 0.0), ga::plane(0.0, 1.0, 0.0, 0.0)), 0, 0, 0, 1, 0, 0,
                   "line_from_planes(x=0, y=0) is the z axis");
    check_vector(ga::plane_through_points(ga::origin<double>(), ga::point(1.0, 0.0, 0.0), ga::point(0.0, 1.0, 0.0)),
                 0, 0, 0, 1, "plane_through_points(O, X, Y) is z = 0");
    check_vector(ga::plane_through_points(ga::origin<double>(), ga::point(1.0, 0.0, 0.0), ga::point(2.0, 0.0, 0.0)),
                 0, 0, 0, 0, "three collinear points give a zero plane");

    // Classification.
    check(!ga::is_ideal(ga::point(1.0, 2.0, 3.0)) && ga::is_ideal(ga::ideal_point(1.0, 2.0, 3.0)), "is_ideal(point)");
    check(!ga::is_ideal(ga::plane(1.0, 0.0, 0.0, 5.0)) && ga::is_ideal(ga::ideal_plane<double>()), "is_ideal(plane)");
    check(!ga::is_ideal(ga::line_through_origin(0.0, 0.0, 1.0)) && ga::is_ideal(Bivector<double>(1, 2, 3, 0, 0, 0)), "is_ideal(line)");
}

void test_geometry()
{
    section("geometry");

    // Reflection in the plane x = 0.
    const Vector<double> mirror = ga::plane(1.0, 0.0, 0.0, 0.0);
    check_close(ga::distance(ga::reflect(ga::point(1.0, 2.0, 3.0), mirror), ga::point(-1.0, 2.0, 3.0)), 0.0,
                "reflecting a point in x = 0 negates its x");
    check(ga::reflect(ga::point(1.0, 2.0, 3.0), mirror).e123 > 0.0, "a reflected point keeps a positive weight");
    // x + y = 1 is plane(1, 1, 0, -1); reflected it becomes -x + y = 1.
    const Vector<double> reflectedPlane = ga::reflect(ga::plane(1.0, 1.0, 0.0, -1.0), mirror);
    check_vector(reflectedPlane, -1, -1, 1, 0, "reflecting x + y = 1 in x = 0 gives -x + y = 1");
    check_bivector(ga::reflect(ga::line_through_origin(0.0, 0.0, 1.0), mirror), 0, 0, 0, 1, 0, 0,
                   "the z axis is unchanged by a reflection in x = 0");
    check_bivector(ga::reflect(ga::line_through_origin(1.0, 0.0, 0.0), mirror), 0, 0, 0, 0, 0, -1,
                   "the x axis is reversed by a reflection in x = 0");
    // Two reflections make a rotation.
    const Vector<double> second = ga::plane(1.0, 1.0, 0.0, 0.0);
    const Trivector<double> twice = ga::reflect(ga::reflect(ga::point(1.0, 0.0, 0.0), mirror), second);
    check_close(ga::distance(twice, ga::sandwich(ga::point(1.0, 0.0, 0.0), ga::motor_between(mirror, second) * ga::motor_between(mirror, second))),
                0.0, "two reflections compose into a rotation by twice the angle");
    // Reflecting in an unnormalised plane is the same as in a normalised one.
    check_close(ga::distance(ga::reflect(ga::point(1.0, 2.0, 3.0), mirror * Scalar<double>(7.0)),
                             ga::reflect(ga::point(1.0, 2.0, 3.0), mirror)), 0.0, "reflect normalises internally");

    // Projection.
    const Vector<double> z1 = ga::plane(0.0, 0.0, 1.0, -1.0);
    check_close(ga::distance(ga::project(ga::point(1.0, 2.0, 3.0), z1), ga::point(1.0, 2.0, 1.0)), 0.0,
                "projecting a point onto z = 1 drops it to z = 1");
    check(ga::project(ga::point(1.0, 2.0, 3.0), z1).e123 > 0.0, "the projected point keeps a positive weight");
    const Bivector<double> zAxis = ga::line_through_origin(0.0, 0.0, 1.0);
    check_close(ga::distance(ga::project(ga::point(1.0, 2.0, 3.0), zAxis), ga::point(0.0, 0.0, 3.0)), 0.0,
                "projecting a point onto the z axis keeps only z");
    check(ga::project(ga::point(1.0, 2.0, 3.0), zAxis).e123 > 0.0, "... with a positive weight");
    const Bivector<double> xLine = ga::line_through_points(ga::point(0.0, -2.0, 0.0), ga::point(1.0, -2.0, 0.0));
    const Bivector<double> projected = ga::project(xLine, z1);
    check_close(ga::distance(ga::point(0.0, -2.0, 1.0), projected), 0.0, "projecting a line onto z = 1 lands it in the plane");
    const std::array<double, 3> projDir = ga::direction(ga::normalize(projected));
    check_close(projDir[0], 1.0, "... keeping its direction, not reversing it");
    // A tilted line must match the line through its two projected points.
    const Trivector<double> diagA = ga::point(0.0, 0.0, 0.0);
    const Trivector<double> diagB = ga::point(1.0, 1.0, 1.0);
    const Vector<double> y3 = ga::plane(0.0, 1.0, 0.0, -3.0);
    check(ga::approx_equal(ga::normalize(ga::project(ga::line_through_points(diagA, diagB), y3)),
                           ga::normalize(ga::line_through_points(ga::project(diagA, y3), ga::project(diagB, y3))),
                           1e-12),
          "projecting a line agrees with joining its projected points, orientation included");
    check_trivector(ga::project(ga::point(1.0, 2.0, 3.0), ga::ideal_plane<double>()), 0, 0, 0, 0,
                    "projecting onto the ideal plane gives zero, not NaN");

    // Distances.
    check_close(ga::distance(ga::point(1.0, 2.0, 3.0), ga::point(4.0, 6.0, 3.0)), 5.0, "distance between two points");
    check_close(ga::distance(ga::point(1.0, 2.0, 3.0), ga::point(1.0, 2.0, 3.0)), 0.0, "a point is zero from itself");
    check_close(ga::distance(ga::point(1.0, 2.0, 3.0) * Scalar<double>(3.0), ga::point(4.0, 6.0, 3.0) * Scalar<double>(-2.0)),
                5.0, "distance divides out the weights");
    check_close(ga::distance(ga::point(1.0, 2.0, 3.0), z1), 2.0, "a point 2 above the plane z = 1");
    check_close(ga::distance(ga::point(1.0, 2.0, -1.0), z1), -2.0, "... and 2 below it, signed");
    check_close(ga::distance(ga::point(1.0, 2.0, 1.0), z1), 0.0, "a point on the plane is zero away");
    check_close(ga::distance(ga::point(1.0, 2.0, 3.0), ga::plane(3.0, 0.0, 0.0, -3.0)), 0.0,
                "distance normalises the plane internally");
    check_close(ga::distance(ga::point(1.0, 2.0, 3.0), zAxis), std::sqrt(5.0), "distance from a point to the z axis");
    check_close(ga::distance(ga::point(0.0, 0.0, 7.0), zAxis), 0.0, "a point on the line is zero away");
    check_close(ga::distance(ga::ideal_point(1.0, 0.0, 0.0), ga::origin<double>()), 0.0, "an ideal point has no distance");

    // Angles.
    check_close(ga::angle(ga::plane(1.0, 0.0, 0.0, 0.0), ga::plane(0.0, 1.0, 0.0, 0.0)), kPi / 2, "perpendicular planes");
    check_close(ga::angle(ga::plane(1.0, 0.0, 0.0, 0.0), ga::plane(1.0, 1.0, 0.0, 0.0)), kPi / 4, "planes at 45 degrees");
    check_close(ga::angle(ga::plane(1.0, 0.0, 0.0, 0.0), ga::plane(1.0, 0.0, 0.0, -9.0)), 0.0,
                "parallel planes are at zero degrees whatever their offsets");
    check_close(ga::angle(zAxis, ga::line_through_origin(1.0, 0.0, 0.0)), kPi / 2, "perpendicular lines");
    check_close(ga::angle(zAxis, ga::line_through_origin(1.0, 0.0, 1.0)), kPi / 4, "lines at 45 degrees");
    check_close(ga::angle(zAxis, ga::line_through_points(ga::point(5.0, 5.0, 0.0), ga::point(5.0, 5.0, 1.0))), 0.0,
                "parallel lines are at zero degrees wherever they are");
    check_close(ga::angle(ga::ideal_plane<double>(), ga::plane(1.0, 0.0, 0.0, 0.0)), 0.0, "an ideal plane has no angle");
}

} // namespace

void run_geometry()
{
    test_primitives();
    test_geometry();
}

} // namespace cctest
