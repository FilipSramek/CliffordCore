// 03 -- Rotations with rotors.
//
// This is what the library is actually good at. Rotors compose by multiplication,
// never gimbal-lock, and interpolate cleanly. There is exactly one sign
// convention that catches everyone, and it gets its own section below.
//
// Build:  ./build_examples.sh 03

#include <cmath>
#include <iostream>

#include "cliffordcore.hpp"
#include "print_helpers.hpp"

using ga::Bivector;
using ga::Multivector;
using ga::Rotor;
using ga::Scalar;
using ga::Vector;

namespace {

constexpr double kPi = 3.14159265358979323846;

const Vector<double> e1(1, 0, 0);
const Vector<double> e2(0, 1, 0);
const Vector<double> e3(0, 0, 1);

void section_basic_rotation()
{
    ex::heading("A quarter turn about +z");

    const Rotor<double> r = ga::rotor_from_axis_angle(e3, kPi / 2);
    ex::print("the rotor", r);
    ex::print("norm(r)", ga::norm(r));
    ex::note("a rotor must be unit length to represent a rotation");

    ex::print("rotate(e1, r)", ga::rotate(e1, r));
    ex::print("rotate(e2, r)", ga::rotate(e2, r));
    ex::print("rotate(e3, r)", ga::rotate(e3, r));
    ex::note("e1 -> e2, e2 -> -e1, and the axis e3 is left alone");
    ex::note("this is the ordinary right-hand rule");

    ex::print("sandwich(e1, r)", ga::sandwich(e1, r));
    ex::note("rotate() is just a friendlier name for sandwich()");
    ex::note("MIND THE ORDER: the vector comes first, the rotor second");

    const Vector<double> v(1, 2, 3);
    ex::print("norm(v)", ga::norm(v));
    ex::print("norm(rotate(v, r))", ga::norm(ga::rotate(v, r)));
    ex::note("rotation preserves length");
}

void section_the_sign_that_bites()
{
    ex::heading("The sign that bites you");

    // Same rotation, two spellings -- and they go opposite ways.
    const Rotor<double> viaExp = ga::exp(Bivector<double>(kPi / 4, 0, 0));
    const Rotor<double> viaAxis = ga::rotor_from_axis_angle(e3, kPi / 2);

    ex::print("exp(pi/4 * e12) on e1", ga::rotate(e1, viaExp));
    ex::print("axis_angle(+z, pi/2) on e1", ga::rotate(e1, viaAxis));
    ex::note("both are quarter turns in the xy plane, and they disagree in sign");
    ex::note("");
    ex::note("three things are going on:");
    ex::note("  1. exp(theta * B) rotates by TWO theta, so pi/4 is a quarter turn");
    ex::note("  2. e12 spans the xy PLANE, so it turns about z -- not about x.");
    ex::note("     pairing xy with the x axis is the classic mistake here");
    ex::note("  3. rotor_from_axis_angle is exp(-(angle/2) * dual(axis)); that");
    ex::note("     minus sign is exactly what makes it the right-hand rule");

    ex::heading("Which plane does which bivector span?");
    ex::print("exp(pi/4 * e12) fixes", ga::rotate(e3, viaExp));
    ex::note("e12 is the xy plane, so e3 is its axis and stays put");

    const Rotor<double> yz = ga::exp(Bivector<double>(0, 0, kPi / 4));
    ex::print("exp(pi/4 * e23) fixes", ga::rotate(e1, yz));
    ex::note("e23 is the yz plane, so e1 is its axis and stays put");
}

void section_composition()
{
    ex::heading("Composing rotations");

    const Rotor<double> aboutZ = ga::rotor_from_axis_angle(e3, kPi / 2);
    const Rotor<double> aboutX = ga::rotor_from_axis_angle(e1, kPi / 2);

    const Rotor<double> twice = aboutZ * aboutZ;
    ex::print("aboutZ * aboutZ", twice);
    ex::note("Rotor * Rotor stays a Rotor -- composition never widens");

    const Vector<double> v(1, 2, 3);
    ex::print("rotate(rotate(v, Z), Z)", ga::rotate(ga::rotate(v, aboutZ), aboutZ));
    ex::print("rotate(v, Z * Z)", ga::rotate(v, twice));
    ex::note("rotating twice == rotating once by the composed rotor");

    ex::print("rotate(e1, Z * X)", ga::rotate(e1, aboutZ * aboutX));
    ex::print("rotate(e1, X * Z)", ga::rotate(e1, aboutX * aboutZ));
    ex::note("order matters: rotations do not commute, and neither does the product");
}

void section_undoing()
{
    ex::heading("Undoing a rotation");

    const Rotor<double> r = ga::rotor_from_axis_angle(Vector<double>(1, 2, 3), 0.9);
    const Vector<double> v(3, -1, 2);

    ex::print("norm(r)", ga::norm(r));
    ex::print("inverse(r)", ga::inverse(r));
    ex::print("reverse(r)", ga::reverse(r));
    ex::note("for a UNIT rotor the inverse is just the reverse -- flip the");
    ex::note("bivector. inverse() divides by the squared norm as well, which");
    ex::note("only matters if the rotor has drifted off unit length");

    ex::print("v", v);
    ex::print("rotate then un-rotate", ga::rotate(ga::rotate(v, r), ga::inverse(r)));
}

void section_the_long_way()
{
    ex::heading("What sandwich() is short for");

    const Rotor<double> r = ga::rotor_from_axis_angle(e3, 1.0);
    const Vector<double> v(1, 2, 3);

    // Rotor * Vector and Multivector * Rotor both give a Multivector.
    const Multivector<double> spelled = r * v * ga::reverse(r);
    ex::multivector_legend();
    ex::print("r * v * reverse(r)", spelled);
    ex::note("every part except the vector one is zero -- the sandwich of a");
    ex::note("vector by a rotor is always another vector");
    ex::print("its vector part", ga::grade1(spelled));
    ex::print("sandwich(v, r)", ga::sandwich(v, r));
    ex::note("identical. sandwich() is the closed form of that product,");
    ex::note("which is why it exists: same answer, far less arithmetic");
}

// Compose `step` with itself `iterations` times and report how far the result
// has drifted from unit length. Templated so we can contrast float and double.
template <typename T>
T accumulated_drift(int iterations)
{
    const ga::Rotor<T> step =
        ga::rotor_from_axis_angle(ga::Vector<T>(0, 0, 1), T(0.01));
    ga::Rotor<T> accumulated = ga::identity_rotor<T>();
    for (int i = 0; i < iterations; ++i) {
        accumulated = accumulated * step;
    }
    return std::abs(ga::norm(accumulated).value - T(1));
}

void section_drift()
{
    ex::heading("Drift, and when normalize() actually matters");

    ex::note("a rotor has to stay unit length; if it drifts it starts scaling");
    ex::note("as well as rotating. Here is |norm - 1| after repeated products:");

    ex::line("double,   1 000 products", ex::sci(accumulated_drift<double>(1000)));
    ex::line("double, 1 000 000 products", ex::sci(accumulated_drift<double>(1000000)));
    ex::note("in double the drift is real but tiny -- a million compositions");
    ex::note("still leaves you accurate to about a part in 10^12");

    ex::line("float,    1 000 products", ex::sci(static_cast<double>(accumulated_drift<float>(1000))));
    ex::line("float,  1 000 000 products", ex::sci(static_cast<double>(accumulated_drift<float>(1000000))));
    ex::note("in float it is a different story: about 1.7% off after a million,");
    ex::note("which is a visible scaling error on every vector you rotate");

    // Show the fix on a rotor that has actually drifted.
    Rotor<float> drifted = ga::identity_rotor<float>();
    const Rotor<float> step = ga::rotor_from_axis_angle(Vector<float>(0, 0, 1), 0.01f);
    for (int i = 0; i < 1000000; ++i) {
        drifted = drifted * step;
    }
    const Vector<float> unit(1, 0, 0);
    ex::line("float: rotated unit length", ex::fmt(ga::norm(ga::rotate(unit, drifted)).value));
    ex::line("  after normalize()", ex::fmt(ga::norm(
                 ga::rotate(unit, ga::normalize(drifted))).value));
    ex::note("so: renormalise periodically when you compose many rotations,");
    ex::note("and especially if you are working in float");
}

} // namespace

int main()
{
    std::cout << "CliffordCore -- 03 rotations\n";

    section_basic_rotation();
    section_the_sign_that_bites();
    section_composition();
    section_undoing();
    section_the_long_way();
    section_drift();

    return 0;
}
