// 04 -- Building rotors from what you actually have, and interpolating them.
//
// You rarely start with a bivector. You start with "point this at that", or
// "halfway between these two orientations". This example covers those, and how
// to get an axis and angle back out again.
//
// Build:  ./build_examples.sh 04

#include <cmath>
#include <iomanip>
#include <iostream>

#include "cliffordcore/cl3.hpp"
#include "print_helpers.hpp"

namespace ga = CliffordCore::Cl3;

using ga::Bivector;
using ga::Rotor;
using ga::Scalar;
using ga::Vector;

namespace {

constexpr double kPi = 3.14159265358979323846;

const Vector<double> e1(1, 0, 0);
const Vector<double> e3(0, 0, 1);

void section_rotor_between()
{
    ex::heading("rotor_between: point one direction at another");

    const Vector<double> from(1, 2, 3);
    const Vector<double> to(-2, 1, 0.5);

    const Rotor<double> r = ga::rotor_between(from, to);
    ex::print("from (normalized)", ga::normalize(from));
    ex::print("to   (normalized)", ga::normalize(to));
    ex::print("rotate(from, r)", ga::rotate(ga::normalize(from), r));
    ex::note("lands exactly on `to` -- that is what the rotor was built for");
    ex::print("norm(r)", ga::norm(r));
    ex::note("inputs need not be unit length; they are normalized internally");
}

void section_degenerate_cases()
{
    ex::heading("The two degenerate cases");

    const Rotor<double> same = ga::rotor_between(e1, Vector<double>(2, 0, 0));
    ex::print("rotor_between(e1, 2*e1)", same);
    ex::print("  applied to e1", ga::rotate(e1, same));
    ex::note("parallel inputs: nothing to do, so the identity rotor");

    const Rotor<double> flip = ga::rotor_between(e1, Vector<double>(-1, 0, 0));
    ex::print("rotor_between(e1, -e1)", flip);
    ex::print("  applied to e1", ga::rotate(e1, flip));
    ex::print("  norm", ga::norm(flip));
    ex::note("antiparallel is the awkward one: EVERY plane containing e1 is a");
    ex::note("valid half turn, so the answer is genuinely ambiguous. The library");
    ex::note("picks a perpendicular plane and gives you a real half turn rather");
    ex::note("than a zero rotor or a NaN");
}

void section_slerp()
{
    ex::heading("slerp: interpolating between orientations");

    const Rotor<double> start = ga::identity_rotor<double>();
    const Rotor<double> end = ga::rotor_from_axis_angle(e3, kPi / 2);

    std::cout << "     t     rotate(e1, slerp)                 angle\n";
    for (int i = 0; i <= 4; ++i) {
        const double t = i * 0.25;
        const Vector<double> v = ga::rotate(e1, ga::slerp(start, end, t));
        const double degrees = std::atan2(v.y, v.x) * 180.0 / kPi;
        std::cout << "   " << std::fixed << std::setprecision(2) << t
                  << "     (" << ex::fmt(v.x) << ", " << ex::fmt(v.y) << ", " << ex::fmt(v.z) << ")"
                  << "     " << std::setw(6) << std::setprecision(1) << degrees << "\n";
    }
    ex::note("evenly spaced in angle, and the endpoints are exact");

    ex::heading("The short way round");
    const Rotor<double> negated(Scalar<double>(-end.scalar.value),
                                 Bivector<double>(-end.bivector.xy, -end.bivector.xz, -end.bivector.yz));
    ex::print("slerp(start, end, 0.5)", ga::rotate(e1, ga::slerp(start, end, 0.5)));
    ex::print("slerp(start, -end, 0.5)", ga::rotate(e1, ga::slerp(start, negated, 0.5)));
    ex::note("a rotor and its negation are the SAME rotation, so slerp flips");
    ex::note("the second operand when they point apart. Without that you would");
    ex::note("get the 270-degree path instead of the 90-degree one");
}

void section_axis_angle_round_trip()
{
    ex::heading("Recovering an axis and angle with log()");

    const Vector<double> axis = ga::normalize(Vector<double>(1, 2, 3));
    const double angle = 0.7;

    const Rotor<double> r = ga::rotor_from_axis_angle(axis, angle);
    const Bivector<double> b = ga::log(r);

    ex::print("original axis", axis);
    ex::print("original angle", Scalar<double>(angle));
    ex::print("log(r)", b);
    ex::print("recovered axis", ga::dual(ga::normalize(b)));
    ex::print("recovered angle", Scalar<double>(2.0 * ga::norm(b).value));
    ex::note("angle = 2*|log(r)| because exp(theta*B) turns by 2*theta,");
    ex::note("and axis = dual(normalize(log(r))): the rotor was built as");
    ex::note("exp(-(angle/2)*dual(axis)), and dual undoes itself up to sign");

    ex::print("exp(log(r))", ga::exp(b));
    ex::print("the original r", r);
    ex::note("exp and log are inverses, so this is a clean round trip");
}

void section_log_edge_cases()
{
    ex::heading("log() on the awkward inputs");

    ex::print("log(identity)", ga::log(ga::identity_rotor<double>()));
    ex::note("no rotation, no plane: the zero bivector");

    // A rotor whose scalar part has drifted a hair outside [-1, 1].
    const Rotor<double> drifted(Scalar<double>(1.0 + 2e-16), Bivector<double>(1, 0, 0));
    const Bivector<double> safe = ga::log(drifted);
    ex::print("log(scalar slightly > 1)", safe);
    ex::note("std::acos is undefined outside [-1,1] and would return NaN here.");
    ex::note("Composition drifts, so this is reachable in practice -- log()");
    ex::note("clamps the scalar before the acos rather than handing back NaN");
}

} // namespace

int main()
{
    std::cout << "CliffordCore -- 04 interpolation\n";

    section_rotor_between();
    section_degenerate_cases();
    section_slerp();
    section_axis_angle_round_trip();
    section_log_edge_cases();

    return 0;
}
