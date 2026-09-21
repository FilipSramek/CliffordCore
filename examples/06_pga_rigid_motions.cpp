// 06 -- Rigid motions in Cl(3,0,1), the projective geometric algebra.
//
// The first five examples are all Cl(3,0), where every object passes through
// the origin and a rotor is the only motion there is. This one switches
// algebra with a single line -- namespace ga = CliffordCore::PGA -- and gets
// planes, lines and points anywhere in space, plus translations and screws,
// out of the same six type names.
//
// What changes with the algebra: a Vector is a plane, a Bivector is a line, a
// Trivector is a point. A fourth basis vector e0 with e0^2 = 0 carries the
// offset from the origin, which is what lets an object sit anywhere.
//
// Build: ./build_examples.sh 06 (or see docs/pga.md)

#include <cmath>
#include <iostream>

#include "cliffordcore/pga.hpp"
#include "print_helpers_pga.hpp"

namespace ga = CliffordCore::PGA;

using ga::Bivector;
using ga::Motor;
using ga::Rotor;
using ga::Translator;
using ga::Trivector;
using ga::Vector;

namespace {

const double kPi = 3.14159265358979323846;

void section_the_three_objects()
{
    ex::heading("Planes, lines and points");

    // Grade 1 is a plane, given by its equation. Grade 3 is a point.
    const Vector<double> ground = ga::plane(0.0, 0.0, 1.0, 0.0);      // z = 0
    const Vector<double> ceiling = ga::plane(0.0, 0.0, 1.0, -3.0);    // z = 3
    const Trivector<double> here = ga::point(1.0, 2.0, 3.0);

    ex::print("ground", ground);
    ex::print("ceiling", ceiling);
    ex::print("point", here);
    ex::note("a Vector is a plane and a Trivector is a point -- the same type");
    ex::note("names as Cl3, different meanings. See docs/pga.md.");

    // Grade 2 is a line, and there are two ways to get one.
    const Bivector<double> where_they_meet = ga::line_from_planes(ground, ga::plane(0.0, 1.0, 0.0, 0.0));
    const Bivector<double> through_two_points = ga::line_through_points(ga::point(1.0, 1.0, 0.0), ga::point(1.0, 1.0, 5.0));

    ex::print("meet(z=0, y=0)", where_they_meet);
    ex::print("join((1,1,0), (1,1,5))", through_two_points);
    ex::note("the meet is the outer product; the join is the regressive one.");
    ex::note("Neither line has to pass through the origin -- that is the point.");

    // Incidence is one product, and it is exactly zero when it should be.
    ex::print("ground ^ point(4,-7,0)", ga::plane(0.0, 0.0, 1.0, 0.0) ^ ga::point(4.0, -7.0, 0.0));
    ex::print("ground ^ point(1,2,3)", ground ^ here);
    ex::note("plane ^ point vanishes exactly when the point lies on the plane.");
}

void section_translation()
{
    ex::heading("Translation, which Cl3 cannot express");

    const Translator<double> east = ga::translator(2.0, 0.0, 0.0);
    ex::print("translator(2,0,0)", east);
    ex::note("1 - (dx e01 + dy e02 + dz e03)/2: the half is what the two-sided");
    ex::note("sandwich needs, the same as a rotor's half angle.");

    // One formula moves every kind of object.
    ex::print("origin moved", ga::sandwich(ga::origin<double>(), east));
    ex::print("plane x=0 moved", ga::sandwich(ga::plane(1.0, 0.0, 0.0, 0.0), east));
    ex::print("z axis moved", ga::sandwich(ga::line_through_origin(0.0, 0.0, 1.0), east));
    ex::note("planes, lines and points all move through M x ~M. In vector");
    ex::note("algebra each of those needs its own special case.");

    // Translators compose by adding their displacements.
    ex::print("east then up", east * ga::translator(0.0, 0.0, 5.0));
    ex::print("... applied to the origin", ga::sandwich(ga::origin<double>(), east * ga::translator(0.0, 0.0, 5.0)));
}

void section_rotation_about_any_line()
{
    ex::heading("Rotation about a line that misses the origin");

    // Cl3's rotation: about an axis through the origin.
    const Rotor<double> about_z = ga::rotor_from_axis_angle(0.0, 0.0, 1.0, kPi / 2);
    ex::print("+90 about +z, applied", ga::sandwich(ga::point(1.0, 0.0, 0.0), about_z));
    ex::note("the same right-hand rule as Cl3: +90 about +z takes (1,0,0) to (0,1,0).");

    // PGA's addition: the axis can be anywhere.
    const Bivector<double> offset_axis = ga::line_through_points(ga::point(1.0, 0.0, 0.0), ga::point(1.0, 0.0, 1.0));
    const Motor<double> half_turn = ga::motor_from_line_angle(offset_axis, kPi);

    ex::print("axis through (1,0,0)", offset_axis);
    ex::print("180 turn: (2,0,0) ->", ga::sandwich(ga::point(2.0, 0.0, 0.0), half_turn));
    ex::print("... and (1,0,5) ->", ga::sandwich(ga::point(1.0, 0.0, 5.0), half_turn));
    ex::note("points on the axis are fixed; (2,0,0) swings round to the origin.");
    ex::note("A Cl3 rotor cannot do this: its axis always passes through 0.");
}

void section_screws()
{
    ex::heading("Screws: every rigid motion is one");

    // Chasles' theorem, as a constructor.
    const Motor<double> rise = ga::screw(ga::line_through_origin(0.0, 0.0, 1.0), kPi / 2, 4.0);
    ex::motor_legend();
    ex::print("screw(z axis, 90, 4)", rise);
    ex::print("(1,0,0) ->", ga::sandwich(ga::point(1.0, 0.0, 0.0), rise));
    ex::note("a quarter turn about z and a rise of 4, in one operation.");

    // A motor recovers its own screw through log().
    const Bivector<double> twist = ga::log(rise);
    ex::print("log(motor)", twist);
    ex::note("log gives back the screw as a bivector -- the axis, with the");
    ex::note("angle as its norm. That bivector is why Twist aliases Bivector.");
    ex::line("exp(log(motor)) == motor", ga::approx_equal(ga::exp(twist), rise) ? "yes" : "no");

    // Composition, and undoing.
    const Motor<double> a = ga::motor_from_line_angle(ga::line_through_origin(1.0, 1.0, 1.0), 0.6);
    const Trivector<double> p = ga::point(2.0, -1.0, 0.5);
    ex::print("(rise * a) applied", ga::sandwich(p, rise * a));
    ex::print("a then rise", ga::sandwich(ga::sandwich(p, a), rise));
    ex::note("motors compose by multiplication, rightmost applied first.");
    ex::print("undone with reverse", ga::sandwich(ga::sandwich(p, rise), ga::reverse(rise)));
}

void section_interpolation()
{
    ex::heading("Interpolating a motion");

    const Motor<double> start = ga::identity_motor<double>();
    const Motor<double> finish = ga::screw(ga::line_through_origin(0.0, 0.0, 1.0), kPi / 2, 4.0);
    const Trivector<double> p = ga::point(1.0, 0.0, 0.0);

    for (int i = 0; i <= 4; ++i) {
        const double t = i / 4.0;
        ex::print("t = " + ex::fmt(t), ga::sandwich(p, ga::slerp(start, finish, t)));
    }
    ex::note("slerp moves along a single screw, so the rotation and the");
    ex::note("translation stay coupled -- the path is a helix, not a corner.");
}

void section_the_classical_answers()
{
    ex::heading("Distance, angle, projection, reflection");

    const Trivector<double> p = ga::point(1.0, 2.0, 3.0);
    const Vector<double> z1 = ga::plane(0.0, 0.0, 1.0, -1.0);
    const Bivector<double> z_axis = ga::line_through_origin(0.0, 0.0, 1.0);

    std::cout << "  point to point            " << ex::fmt(ga::distance(p, ga::point(4.0, 6.0, 3.0))) << "\n";
    std::cout << "  point to plane z=1        " << ex::fmt(ga::distance(p, z1)) << "\n";
    std::cout << "  point to the z axis       " << ex::fmt(ga::distance(p, z_axis)) << "\n";
    ex::note("the point-plane distance is signed: which side you are on.");

    std::cout << "  angle between two planes  " << ex::fmt(ga::angle(ga::plane(1.0, 0.0, 0.0, 0.0), ga::plane(1.0, 1.0, 0.0, -8.0))) << " rad\n";
    ex::note("the offsets -- 0 and -8 -- do not affect the angle, because e0");
    ex::note("squares to zero and drops out of the inner product.");

    ex::print("project onto z = 1", ga::project(p, z1));
    ex::print("project onto the z axis", ga::project(p, z_axis));
    ex::print("reflect in x = 0", ga::reflect(p, ga::plane(1.0, 0.0, 0.0, 0.0)));
    ex::print("reflect in x = 3", ga::reflect(p, ga::plane(1.0, 0.0, 0.0, -3.0)));
    ex::note("the two mirrors are parallel -- same normal, different offset --");
    ex::note("and give different answers. A normal vector alone cannot tell");
    ex::note("them apart; the e0 component is what does.");
}

} // namespace

int main()
{
    std::cout << "CliffordCore example 06 -- rigid motions in Cl(3,0,1)\n";

    section_the_three_objects();
    section_translation();
    section_rotation_about_any_line();
    section_screws();
    section_interpolation();
    section_the_classical_answers();

    std::cout << "\n";
    return 0;
}
