// Cl(3,0) tests -- comparison_and_printing: exact == and != and approx_equal.
// Cl3's to_string checks live in instantiation_sweep.cpp, alongside the call
// that instantiates every type's printing.
//
// One category of the Cl3 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

void test_comparison()
{
    section("comparison");

    const Vector<double> a(1, 2, 3);
    const Vector<double> b(1, 2, 3);
    const Vector<double> c(1, 2, 4);

    check(a == b, "identical vectors compare equal");
    check(!(a == c), "differing vectors do not");
    check(a != c, "operator!= is the negation");
    check(!(a != b), "operator!= agrees with operator==");

    check(Scalar<double>(2) == Scalar<double>(2), "Scalar equality");
    check(Bivector<double>(1, 2, 3) == Bivector<double>(1, 2, 3), "Bivector equality");
    check(Trivector<double>(4) == Trivector<double>(4), "Trivector equality");
    check(make_mv(1, 2, 3, 4, 5, 6, 7, 8) == make_mv(1, 2, 3, 4, 5, 6, 7, 8), "Multivector equality");
    check(!(make_mv(1, 2, 3, 4, 5, 6, 7, 8) == make_mv(1, 2, 3, 4, 5, 6, 7, 9)),
          "Multivector inequality catches one differing component");

    // Exact equality is too strict once trigonometry is involved -- which is
    // exactly why approx_equal is a separate function rather than a fuzzy ==.
    const Rotor<double> r = ga::rotor_from_axis_angle(Vector<double>(0, 0, 1), kPi / 2);
    const Vector<double> turned = ga::rotate(Vector<double>(1, 0, 0), r);
    const Vector<double> expected(0, 1, 0);
    check(!(turned == expected), "a rotation leaves residue, so exact == fails");
    check(ga::approx_equal(turned, expected), "approx_equal accepts it");
    check(ga::approx_equal(turned, expected, 1e-9), "approx_equal takes an explicit tolerance");
    check(!ga::approx_equal(turned, expected, 1e-30), "and honours a strict one");

    // A rotor and its negation are the same rotation but different values.
    const Rotor<double> negated(Scalar<double>(-r.scalar.value),
                                 Bivector<double>(-r.bivector.xy, -r.bivector.xz, -r.bivector.yz));
    check(!(r == negated), "a rotor does not equal its negation");
    check(ga::approx_equal(ga::rotate(Vector<double>(1, 0, 0), r),
                                     ga::rotate(Vector<double>(1, 0, 0), negated)),
          "...even though both describe the same rotation");
}

// Exact equality stays usable in a constant expression; a tolerant one could not.
constexpr Vector<double> kCmpA(1, 0, 0);
constexpr Vector<double> kCmpB(1, 0, 0);
static_assert(kCmpA == kCmpB, "operator== is constexpr");
static_assert(!(kCmpA != kCmpB), "operator!= is constexpr");

} // namespace

void run_comparison_and_printing()
{
    test_comparison();
}

} // namespace cctest
