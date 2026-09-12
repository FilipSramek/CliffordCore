// 01 -- Getting started: the three products.
//
// Geometric algebra gives you one product that subsumes the dot and cross
// products you already know. This example builds it up from the two halves.
//
// Build:  ./build_examples.sh 01      (or see docs/getting-started.md)

#include <iostream>

#include "cliffordcore/cl3.hpp"
#include "print_helpers.hpp"

namespace ga = CliffordCore::Cl3;

using ga::Bivector;
using ga::Multivector;
using ga::Rotor;
using ga::Scalar;
using ga::Vector;

namespace {

void section_vectors()
{
    ex::heading("Vectors");

    const Vector<double> a(1, 2, 3);
    const Vector<double> b(4, 5, 6);

    ex::print("a", a);
    ex::print("b", b);
    ex::print("a + b", a + b);
    ex::print("-a", -a);
}

void section_dot()
{
    ex::heading("Dot product  (grade 1 x grade 1 -> grade 0)");

    const Vector<double> a(1, 2, 3);
    const Vector<double> b(4, 5, 6);

    ex::print("a | b", a | b);
    ex::print("dot_product(a, b)", ga::dot_product(a, b));
    ex::note("the operator and the named function are the same thing");

    const Vector<double> e1(1, 0, 0);
    const Vector<double> e2(0, 1, 0);
    ex::print("e1 | e2", e1 | e2);
    ex::note("zero: perpendicular vectors have no shared direction");
}

void section_wedge()
{
    ex::heading("Wedge product  (grade 1 x grade 1 -> grade 2)");

    const Vector<double> a(1, 2, 3);
    const Vector<double> b(4, 5, 6);

    ex::print("a ^ b", a ^ b);
    ex::note("a bivector: the oriented plane the two vectors span");
    ex::note("components are (xy, xz, yz) = (e1e2, e1e3, e2e3) -- note xz, not zx");

    ex::print("b ^ a", b ^ a);
    ex::note("reversing the operands flips the orientation: b^a == -(a^b)");

    ex::print("a ^ a", a ^ a);
    ex::note("zero: a vector spans no plane with itself");
}

void section_geometric()
{
    ex::heading("Geometric product  (the two halves together)");

    const Vector<double> a(1, 2, 3);
    const Vector<double> b(4, 5, 6);

    const Multivector<double> product = a * b;
    ex::multivector_legend();
    ex::print("a * b", product);
    ex::print("  its scalar part", product.scalar);
    ex::print("  its bivector part", product.bivector);
    ex::note("a*b = (a|b) + (a^b) -- the dot and the wedge, side by side");
    ex::note("that is the whole idea: one product keeping both pieces");

    ex::print("a * a", a * a);
    ex::note("a vector times itself is its squared length, and nothing else");
}

void section_return_types()
{
    ex::heading("Return types");

    const Vector<double> a(1, 2, 3);
    const Vector<double> b(4, 5, 6);

    const Multivector<double> wide = ga::geometric_product(a, b);
    const Rotor<double> narrow = ga::rotor_product(a, b);

    ex::print("geometric_product(a, b)", wide);
    ex::print("rotor_product(a, b)", narrow);
    ex::note("identical numbers, narrower type -- a Rotor holds exactly the");
    ex::note("scalar and bivector parts, which is all this product ever produces");
    ex::note("operator* cannot offer both: C++ cannot overload on return type,");
    ex::note("so operator* gives the general type and the narrow one is named");
}

void section_magnitudes()
{
    ex::heading("Magnitudes");

    const Vector<double> a(1, 2, 3);

    ex::print("norm(a)", ga::norm(a));
    ex::print("squared_norm(a)", ga::squared_norm(a));
    ex::note("squared_norm avoids the square root when you only need a comparison");

    const Vector<double> unit = ga::normalize(a);
    ex::print("normalize(a)", unit);
    ex::print("norm(normalize(a))", ga::norm(unit));
}

void section_scaling()
{
    ex::heading("Scaling -- three spellings, one result");

    const Vector<double> a(1, 2, 3);

    ex::print("a * Scalar<double>(2)", a * Scalar<double>(2));
    ex::print("a * 2.0", a * 2.0);
    ex::print("2.0 * a", 2.0 * a);

    Vector<double> mutable_copy = a;
    mutable_copy *= 2.0;
    ex::print("a *= 2.0", mutable_copy);
}

// Everything above also works at compile time. If any of these were wrong the
// program would not build, so they cost nothing at runtime.
constexpr Vector<double> kE1(1, 0, 0);
constexpr Vector<double> kE2(0, 1, 0);
static_assert((kE1 | kE2).value == 0.0, "perpendicular vectors have zero dot product");
static_assert((kE1 ^ kE2).xy == 1.0, "e1 ^ e2 is the unit xy bivector");
static_assert((kE1 * kE2).bivector.xy == 1.0, "the geometric product agrees");
static_assert((kE1 * kE1).scalar.value == 1.0, "a unit vector squares to 1");

} // namespace

int main()
{
    std::cout << "CliffordCore -- 01 getting started\n";

    section_vectors();
    section_dot();
    section_wedge();
    section_geometric();
    section_return_types();
    section_magnitudes();
    section_scaling();

    ex::heading("Compile time");
    ex::note("the static_asserts in this file proved the products above at");
    ex::note("compile time -- the whole library is constexpr");

    return 0;
}
