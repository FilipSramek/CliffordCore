// 05 -- The classical 3D geometry answers, written in GA.
//
// Angles, areas, volumes, projections and reflections. Each one is a one-liner
// once you have the products, and the last section shows why rotors exist at all.
//
// Build:  ./build_examples.sh 05

#include <cmath>
#include <iostream>

#include "cliffordcore.hpp"
#include "print_helpers.hpp"

using CliffordCore::Bivector3;
using CliffordCore::Multivector3;
using CliffordCore::Rotor3;
using CliffordCore::Scalar;
using CliffordCore::Trivector3;
using CliffordCore::Vector3;

namespace {

constexpr double kPi = 3.14159265358979323846;

void section_angle()
{
    ex::heading("Angle between two vectors");

    const Vector3<double> a(1, 0, 0);
    const Vector3<double> b(1, 1, 0);

    // The dot gives |a||b|cos, the wedge's norm gives |a||b|sin. Feeding both
    // to atan2 is better conditioned than acos of the normalized dot, which
    // loses precision badly for nearly-parallel vectors.
    const double radians = std::atan2(CliffordCore::norm(a ^ b).value, (a | b).value);

    ex::print("a", a);
    ex::print("b", b);
    ex::print("angle (degrees)", Scalar<double>(radians * 180.0 / kPi));
    ex::note("atan2(|a^b|, a|b): the wedge supplies the sine, the dot the cosine");
    ex::note("better conditioned than acos() when the vectors are near-parallel");
}

void section_area_and_volume()
{
    ex::heading("Area and volume");

    const Vector3<double> a(3, 0, 0);
    const Vector3<double> b(0, 4, 0);
    const Vector3<double> c(0, 0, 5);

    ex::print("|a ^ b| (parallelogram)", CliffordCore::norm(a ^ b));
    ex::print("triangle area", Scalar<double>(CliffordCore::norm(a ^ b).value / 2.0));

    const Trivector3<double> volume = (a ^ b) ^ c;
    ex::print("(a ^ b) ^ c", volume);
    ex::note("the signed volume of the box, i.e. the scalar triple product");
    ex::note("this is the only geometric route to a trivector in the library:");
    ex::note("wedging a plane with a vector outside it");

    ex::print("(a ^ b) ^ a", (a ^ b) ^ a);
    ex::note("zero: a lies in the plane a^b, so the box is flat");
}

void section_projection()
{
    ex::heading("Projection and rejection");

    const Vector3<double> v(2, 3, 4);
    const Vector3<double> u = CliffordCore::normalize(Vector3<double>(1, 1, 0));

    const Vector3<double> parallel = u * (v | u);
    const Vector3<double> perpendicular = v - parallel;

    ex::print("v", v);
    ex::print("u (unit)", u);
    ex::print("projection onto u", parallel);
    ex::print("rejection from u", perpendicular);
    ex::print("  the two summed", parallel + perpendicular);
    ex::print("  their dot (want 0)", parallel | perpendicular);

    // The same rejection, spelled with the wedge: the part of v that spans a
    // plane with u, divided back through u.
    const Vector3<double> viaWedge = CliffordCore::grade1((v ^ u) * CliffordCore::inverse(u));
    ex::print("rejection via (v^u)/u", viaWedge);
    ex::note("the wedge form generalises: it works for any invertible u,");
    ex::note("and to planes as well as vectors");
}

void section_reflection()
{
    ex::heading("Reflection in a plane");

    const Vector3<double> n = CliffordCore::normalize(Vector3<double>(0, 1, 0));
    const Vector3<double> v(2, 3, 4);

    // In GA, reflecting v in the plane with unit normal n is -n v n.
    const Vector3<double> reflected = CliffordCore::grade1(-(n * v * n));
    const Vector3<double> classical = v - 2.0 * (n * (v | n));

    ex::print("v", v);
    ex::print("plane normal n", n);
    ex::print("-n v n", reflected);
    ex::print("v - 2(v.n)n", classical);
    ex::note("the same reflection, one written as a product and one by hand");
}

void section_two_reflections()
{
    ex::heading("Two reflections make a rotation");

    // Two planes through the origin, 45 degrees apart in the xy plane.
    const Vector3<double> n = CliffordCore::normalize(Vector3<double>(1, 0, 0));
    const Vector3<double> m = CliffordCore::normalize(Vector3<double>(1, 1, 0));
    const Vector3<double> v(1, 2, 3);

    const Vector3<double> once = CliffordCore::grade1(-(n * v * n));
    const Vector3<double> twice = CliffordCore::grade1(-(m * once * m));

    // ...and the same thing as a single rotor.
    const Rotor3<double> r = CliffordCore::rotor_product(m, n);
    const Vector3<double> viaRotor = CliffordCore::rotate(v, r);

    ex::print("v", v);
    ex::print("reflect in n, then in m", twice);
    ex::print("rotate(v, rotor_product(m,n))", viaRotor);
    ex::note("identical. Reflecting in two planes IS a rotation -- by twice");
    ex::note("the angle between them, about their line of intersection");
    ex::note("");
    ex::note("that is what a rotor is: the product of two unit vectors, which");
    ex::note("is why rotor_product(a, b) exists and why the sandwich R v ~R");
    ex::note("has the shape it does -- it is two reflections in a row");

    const double planeAngle = std::atan2(CliffordCore::norm(m ^ n).value, (m | n).value);
    ex::print("angle between planes (deg)", Scalar<double>(planeAngle * 180.0 / kPi));
    ex::print("rotation angle (deg)", Scalar<double>(2.0 * planeAngle * 180.0 / kPi));
}

void section_inverse_check()
{
    ex::heading("A vector times its inverse");

    const Vector3<double> v(1, 2, 3);
    ex::multivector_legend();
    ex::print("v * inverse(v)", v * CliffordCore::inverse(v));
    ex::note("exactly 1, nothing else -- as it should be");
}

constexpr Vector3<double> kE1(1, 0, 0);
constexpr Vector3<double> kE2(0, 1, 0);
constexpr Vector3<double> kE3(0, 0, 1);
static_assert(((kE1 ^ kE2) ^ kE3).e123 == 1.0, "the unit box has volume +1");
static_assert(((kE1 ^ kE2) ^ kE1).e123 == 0.0, "a flat box has zero volume");

} // namespace

int main()
{
    std::cout << "CliffordCore -- 05 geometry toolkit\n";

    section_angle();
    section_area_and_volume();
    section_projection();
    section_reflection();
    section_two_reflections();
    section_inverse_check();

    return 0;
}
