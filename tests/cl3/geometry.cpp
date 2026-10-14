// Cl(3,0) tests -- geometry: reflect, project and reject.
//
// One category of the Cl3 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

void test_geometry()
{
    section("geometry");

    const Vector<double> v(2, 3, 4);
    const Vector<double> n(0, 1, 0);
    const Bivector<double> xy(1, 0, 0);

    // Reflection matches the formula everyone already knows.
    check_vector(ga::reflect(v, n), 2.0, -3.0, 4.0, "reflect in the plane normal to e2");
    const Vector<double> classical = v - 2.0 * (n * (v | n));
    check_vector(ga::reflect(v, n), classical.x, classical.y, classical.z,
                 "reflect matches v - 2(v.n)n");

    // The normal need not be unit length -- reflect divides through by inverse(n).
    check_vector(ga::reflect(v, Vector<double>(0, 5, 0)), 2.0, -3.0, 4.0,
                 "reflect normalizes the plane normal itself");

    // Reflecting twice is the identity.
    check_vector(ga::reflect(ga::reflect(v, n), n), v.x, v.y, v.z,
                 "reflect twice returns the original");

    // A plane can be given as a bivector rather than by its normal.
    check_vector(ga::reflect(v, xy), 2.0, 3.0, -4.0, "reflect in the e12 plane");

    // Projection and rejection split a vector in two, onto a line...
    const Vector<double> u(1, 1, 0);
    const Vector<double> par = ga::project(v, u);
    const Vector<double> perp = ga::reject(v, u);
    check_vector(par, 2.5, 2.5, 0.0, "project onto a line");
    check_vector(par + perp, v.x, v.y, v.z, "project + reject == v, for a line");
    check_scalar(perp | u, 0.0, "the rejection is perpendicular to the line");

    // ...and onto a plane, which is what needed the contraction.
    const Vector<double> inPlane = ga::project(v, xy);
    const Vector<double> outOfPlane = ga::reject(v, xy);
    check_vector(inPlane, 2.0, 3.0, 0.0, "project onto the e12 plane drops the z part");
    check_vector(outOfPlane, 0.0, 0.0, 4.0, "reject from the e12 plane keeps only z");
    check_vector(inPlane + outOfPlane, v.x, v.y, v.z, "project + reject == v, for a plane");

    // The same, on a plane that is not axis aligned.
    const Bivector<double> tilted(1, 2, 3);
    const Vector<double> tiltedSum = ga::project(v, tilted) + ga::reject(v, tilted);
    check_vector(tiltedSum, v.x, v.y, v.z, "project + reject == v, tilted plane", 1e-14);
    const Vector<double> reprojected = ga::project(ga::reject(v, tilted), tilted);
    check_vector(reprojected, 0.0, 0.0, 0.0, "the rejection has nothing left in the plane", 1e-14);
}

} // namespace

void run_geometry()
{
    test_geometry();
}

} // namespace cctest
