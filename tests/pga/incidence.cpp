// Cl(3,0,1) tests -- incidence: incidence, join and meet: the outer and
// regressive products on points, lines and planes.
//
// One category of the PGA suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

void test_incidence()
{
    section("incidence");

    // plane ^ point = (a x + b y + c z + d) e0123, verified in docs/pga.md.
    const Vector<double> p = test_plane(1, 2, 3, 4);
    const Trivector<double> P = test_point(0.5, -1, 2);
    check_quadvector(p ^ P, 1 * 0.5 + 2 * -1 + 3 * 2 + 4, "plane ^ point is the plane equation");
    check_quadvector(P ^ p, -(1 * 0.5 + 2 * -1 + 3 * 2 + 4), "point ^ plane has the opposite sign");

    // A point on the plane gives exactly zero.
    const Vector<double> x1 = test_plane(1, 0, 0, -1);
    check((x1 ^ test_point(1, 7, -3)).e0123 == 0.0, "point(1, 7, -3) lies on x = 1");
    check((x1 ^ test_point(2, 0, 0)).e0123 != 0.0, "point(2, 0, 0) does not lie on x = 1");

    // The origin is e123 and lies on every plane through the origin.
    check((test_plane(3, -2, 5, 0) ^ test_point(0, 0, 0)).e0123 == 0.0, "the origin lies on 3x - 2y + 5z = 0");

    // The point is the triple meet of its coordinate planes.
    const Vector<double> px = test_plane(1, 0, 0, -1), py = test_plane(0, 1, 0, -2), pz = test_plane(0, 0, 1, -3);
    check_trivector((px ^ py) ^ pz, -3, 2, -1, 1, "point(1,2,3) == (x=1) ^ (y=2) ^ (z=3)");
    check_trivector(px ^ (py ^ pz), -3, 2, -1, 1, "the triple meet associates");

    // Two lines: zero iff coplanar.
    const Bivector<double> zAxis(0, 0, 0, 1, 0, 0);
    const Bivector<double> zParallel(0, -1, 0, 1, 0, 0);   // through (1,0,0)
    const Bivector<double> xAxis(0, 0, 0, 0, 0, 1);
    check((zAxis ^ zParallel).e0123 == 0.0, "parallel lines are coplanar");
    check((zAxis ^ xAxis).e0123 == 0.0, "intersecting lines are coplanar");
    const Bivector<double> skew = test_plane(0, 1, 0, -1) ^ test_plane(0, 0, 1, 0);   // y = 1, z = 0: parallel to x, offset
    check((zAxis ^ skew).e0123 != 0.0, "skew lines have a non-zero wedge");
}

void test_join_and_meet()
{
    section("join and meet");

    const Trivector<double> origin = test_point(0, 0, 0);

    // The lines through the origin, from docs/pga.md.
    check_bivector(ga::join(origin, test_point(1, 0, 0)), 0, 0, 0, 0, 0, 1, "join(origin, x) == e23");
    check_bivector(ga::join(origin, test_point(0, 1, 0)), 0, 0, 0, 0, -1, 0, "join(origin, y) == -e13");
    check_bivector(ga::join(origin, test_point(0, 0, 1)), 0, 0, 0, 1, 0, 0, "join(origin, z) == e12");
    check_bivector(ga::join(test_point(1, 0, 0), test_point(1, 0, 1)), 0, -1, 0, 1, 0, 0,
                   "the z-parallel line through (1,0,0) is e12 - e02");
    check_bivector(ga::join(test_point(1, 0, 1), test_point(1, 0, 0)), 0, 1, 0, -1, 0, 0, "join is antisymmetric on points");
    check_bivector(ga::join(origin, origin), 0, 0, 0, 0, 0, 0, "join of a point with itself is zero");

    // Distance between points is the Euclidean norm of their join.
    check_scalar(ga::norm(ga::join(test_point(1, 2, 3), test_point(4, 6, 3))), 5.0, "norm(join(P, Q)) is the distance");

    // Point-line join is the plane through both; zero iff incident.
    const Bivector<double> zAxis(0, 0, 0, 1, 0, 0);
    const Vector<double> through = ga::join(test_point(1, 2, 3), zAxis);
    check_scalar(ga::norm(through), std::sqrt(5.0), "norm(join(P, L)) is the distance from P to L");
    check((through ^ test_point(1, 2, 3)).e0123 == 0.0, "the joined plane contains the point");
    check((through ^ origin).e0123 == 0.0, "... and the line");
    check_vector(ga::join(test_point(0, 0, 5), zAxis), 0, 0, 0, 0, "a point on the line joins to zero");
    // dual(L) ^ dual(P) is a bivector wedge a vector, which commutes.
    check_vector(ga::join(zAxis, test_point(1, 2, 3)), through.e0, through.e1, through.e2, through.e3,
                 "join(line, point) == join(point, line)");

    // Plane through three points.
    const Vector<double> z0 = ga::join(ga::join(origin, test_point(1, 0, 0)), test_point(0, 1, 0));
    check_vector(z0, 0, 0, 0, 1, "join(join(O, X), Y) is the plane z = 0, e3");
    check((z0 ^ test_point(4, -7, 0)).e0123 == 0.0, "the plane through O, X, Y contains (4, -7, 0)");

    // Two lines: zero iff coplanar.
    check_scalar(ga::join(zAxis, Bivector<double>(0, -1, 0, 1, 0, 0)), 0.0, "parallel lines join to zero");
    check_scalar(ga::join(zAxis, Bivector<double>(0, 0, 0, 0, 0, 1)), 0.0, "intersecting lines join to zero");
    const Bivector<double> skew = test_plane(0, 1, 0, -1) ^ test_plane(0, 0, 1, 0);
    check(ga::join(zAxis, skew).value != 0.0, "skew lines join to a non-zero scalar");

    // meet is the wedge under its geometric name.
    check_bivector(ga::meet(test_plane(1, 0, 0, 0), test_plane(0, 1, 0, 0)), 0, 0, 0, 1, 0, 0, "meet(x=0, y=0) is the z-axis");
    const Trivector<double> hit = ga::meet(test_plane(0, 0, 1, -1), zAxis);
    check_trivector(ga::normalize(hit), -1, 0, 0, 1, "meet(z=1, z-axis) is point(0, 0, 1)");
    check_trivector(ga::meet(zAxis, test_plane(0, 0, 1, -1)), hit.e012, hit.e013, hit.e023, hit.e123, "meet(line, plane) == meet(plane, line)");

    // The named forms are the same thing.
    check_bivector(ga::regressive_product(origin, test_point(0, 0, 1)), 0, 0, 0, 1, 0, 0, "regressive_product == join");
    check_scalar(ga::regressive_product(zAxis, skew), ga::join(zAxis, skew).value, "regressive_product(line, line) == join");
}

} // namespace

void run_incidence()
{
    test_incidence();
    test_join_and_meet();
}

} // namespace cctest
