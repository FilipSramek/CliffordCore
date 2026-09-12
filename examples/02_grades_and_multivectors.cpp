// 02 -- Grades, multivectors, and duality.
//
// Cl(3,0) has four grades: scalars, vectors, bivectors (planes) and the
// pseudoscalar (volume). This example shows how they combine, how to pull them
// apart again, and the two places the type system will surprise you.
//
// Build:  ./build_examples.sh 02

#include <iostream>

#include "cliffordcore/cl3.hpp"
#include "print_helpers.hpp"

namespace ga = CliffordCore::Cl3;

using ga::Bivector;
using ga::Multivector;
using ga::Rotor;
using ga::Scalar;
using ga::Trivector;
using ga::Vector;

namespace {

void section_the_four_grades()
{
    ex::heading("The four grades");

    ex::print("grade 0  Scalar", Scalar<double>(2));
    ex::print("grade 1  Vector", Vector<double>(1, 2, 3));
    ex::print("grade 2  Bivector", Bivector<double>(1, 2, 3));
    ex::print("grade 3  Trivector", Trivector<double>(4));
    ex::note("grade 3 is the pseudoscalar: the oriented volume element");
    ex::note("its member is .e123 -- NOT .value, which is Scalar's member");
}

void section_widening()
{
    ex::heading("Adding across grades widens to a Multivector");

    const Scalar<double> s(2);
    const Vector<double> v(1, 2, 3);
    const Bivector<double> b(4, 5, 6);
    const Trivector<double> t(7);

    ex::multivector_legend();
    ex::print("s + v", s + v);
    ex::print("v + b", v + b);
    ex::print("b + t", b + t);
    ex::note("different grades cannot share a type, so the sum widens");

    ex::print("s + v + b + t", s + v + b + t);
    ex::note("chains work because a+b+c parses as ((a+b)+c) -- only pairs");
    ex::note("are defined, and every chain length follows from that");

    ex::print("v + v", v + v);
    ex::note("...but a same-grade sum stays in its own grade: still a Vector");

    ex::heading("The one case that looks like it should be special");
    ex::print("s + b", s + b);
    ex::note("a scalar plus a bivector IS exactly a rotor, yet this still");
    ex::note("widens -- uniformity keeps chaining predictable, and C++ cannot");
    ex::note("overload on return type. Ask for a rotor by name instead:");
    ex::print("rotor_sum(s, b)", ga::rotor_sum(s, b));
    ex::print("to_rotor(s + b)", ga::to_rotor(s + b));
    ex::print("to_multivector(rotor)", ga::to_multivector(ga::rotor_sum(s, b)));
}

void section_projection()
{
    ex::heading("Pulling the grades back out");

    const Multivector<double> m(Scalar<double>(1), Vector<double>(2, 3, 4),
                                 Bivector<double>(5, 6, 7), Trivector<double>(8));
    ex::multivector_legend();
    ex::print("m", m);
    ex::print("grade0(m)", ga::grade0(m));
    ex::print("grade1(m)", ga::grade1(m));
    ex::print("grade2(m)", ga::grade2(m));
    ex::print("grade3(m)", ga::grade3(m));
    ex::note("named per grade rather than grade<N>(m) because each returns a");
    ex::note("different type, which one template cannot express cleanly");

    const Multivector<double> rebuilt = ga::grade0(m) + ga::grade1(m)
                                       + ga::grade2(m) + ga::grade3(m);
    ex::print("the four, added back", rebuilt);
}

void section_trivector_trap()
{
    ex::heading("The Trivector return-type trap");

    const Trivector<double> t(6);
    const Trivector<double> u(3);

    ex::print("t", t);
    ex::print("u", u);
    ex::print("t * u   -> Scalar", t * u);
    ex::note("two pseudoscalars multiply to a SCALAR, not a trivector,");
    ex::note("because e123 * e123 = -1. So 6 * 3 comes out as -18");
    ex::print("t / u   -> Scalar", t / u);
    ex::note("division is t * inverse(u), and the two minus signs cancel: 6/3");

    ex::print("t * 2.0 -> Trivector", t * 2.0);
    ex::print("t * Scalar(2) -> Trivector", t * Scalar<double>(2));
    ex::note("but SCALING a pseudoscalar leaves it a pseudoscalar");
    ex::note("so t*u and t*2.0 deliberately differ in return type");
}

void section_duality()
{
    ex::heading("Duality: multiplying by the pseudoscalar");

    const Vector<double> v(1, 2, 3);

    ex::print("dual(Scalar(7))", ga::dual(Scalar<double>(7)));
    ex::print("dual(v)", ga::dual(v));
    ex::print("dual(Bivector(4,5,6))", ga::dual(Bivector<double>(4, 5, 6)));
    ex::print("dual(Trivector(7))", ga::dual(Trivector<double>(7)));
    ex::note("the dual swaps grade k with grade 3-k");

    ex::print("dual(dual(v))", ga::dual(ga::dual(v)));
    ex::note("applying it twice NEGATES rather than returning v -- that is");
    ex::note("correct, not a bug: the pseudoscalar squares to -1");

    const Vector<double> a(1, 0, 0);
    const Vector<double> b(0, 1, 0);
    ex::print("dual(b ^ a)", ga::dual(b ^ a));
    ex::note("this is the familiar cross product a x b -- the plane a^b and");
    ex::note("the axis a x b are dual descriptions of the same thing");
}

void section_involutions()
{
    ex::heading("Reverse and conjugate");

    const Multivector<double> m(Scalar<double>(1), Vector<double>(2, 3, 4),
                                 Bivector<double>(5, 6, 7), Trivector<double>(8));
    ex::multivector_legend();
    ex::print("m", m);
    ex::print("reverse(m)", ga::reverse(m));
    ex::note("grades 0 and 1 keep their sign, grades 2 and 3 flip");
    ex::print("conjugate(m)", ga::conjugate(m));
    ex::note("Clifford conjugation instead flips grades 1 and 2");
}

void section_inverses()
{
    ex::heading("Inverses");

    const Vector<double> v(1, 2, 3);
    const Bivector<double> b(1, 2, 3);
    const Trivector<double> t(4);

    ex::print("inverse(v)", ga::inverse(v));
    ex::print("inverse(b)", ga::inverse(b));
    ex::note("note the minus sign: a bivector squares to -|b|^2");
    ex::print("inverse(t)", ga::inverse(t));
    ex::note("same reason -- the pseudoscalar also squares to -1");

    const Multivector<double> m(Scalar<double>(1), Vector<double>(2, -3, 4),
                                 Bivector<double>(-5, 6, 7), Trivector<double>(-8));
    ex::multivector_legend();
    ex::print("m", m);
    ex::print("inverse(m)", ga::inverse(m));
    ex::print("m * inverse(m)", m * ga::inverse(m));
    ex::note("exactly 1, with every other grade zero");
    ex::note("this is NOT m/|m|^2 -- that only works for a single grade.");
    ex::note("m * conjugate(m) collapses to a scalar plus a pseudoscalar,");
    ex::note("which inverts like a complex number since e123^2 = -1");
}

constexpr Vector<double> kA(1, 0, 0);
constexpr Vector<double> kB(0, 1, 0);
static_assert(ga::dual(kB ^ kA).z == 1.0, "dual(b^a) is the cross product a x b");
static_assert(ga::dual(ga::dual(kA)).x == -1.0, "dual twice negates");
// Compare .value explicitly: the library has no operator== (see the known gaps
// in docs/conventions.md), so there is nothing to compare Scalar against.
static_assert((Trivector<double>(1) * Trivector<double>(1)).value == -1.0,
              "e123 * e123 == -1");

} // namespace

int main()
{
    std::cout << "CliffordCore -- 02 grades and multivectors\n";

    section_the_four_grades();
    section_widening();
    section_projection();
    section_trivector_trap();
    section_duality();
    section_involutions();
    section_inverses();

    return 0;
}
