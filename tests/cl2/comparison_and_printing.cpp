// Cl(2,0) tests -- comparison_and_printing: to_string, exact == and !=,
// approx_equal, and operator<<.
//
// One category of the Cl2 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

// ---------------------------------------------------------------------------
// String formatting
// ---------------------------------------------------------------------------

void test_to_string()
{
    section("to_string");

    check(Scalar<double>(1.5).to_string() == "1.5", "Scalar::to_string");
    check(Vector<double>(1.0, 2.0).to_string() == "1*e1 + 2*e2", "Vector::to_string");
    check(Bivector<double>(3.0).to_string() == "3*e12", "Bivector::to_string");
    check(Rotor<double>(Scalar<double>(1.0), Bivector<double>(2.0)).to_string()
              == "(1) + (2*e12)",
          "Rotor::to_string");
    check(make_mv(1.0, 2.0, 3.0, 4.0).to_string() == "(1) + (2*e1 + 3*e2) + (4*e12)",
          "Multivector::to_string");

    // Round-trip precision, not std::to_string's fixed six decimals: a 1e-16
    // residual must not print as 0.000000.
    check(Scalar<double>(1e-16).to_string() != "0.000000", "to_string keeps small values");
}

void test_comparison()
{
    section("comparison");

    // The default tolerance is 100 epsilons, about 2.2e-14 for double: enough
    // to absorb rounding, and deliberately too small to hide a 1e-10 error.
    const Vector<double> a(1.0, 2.0), b(1.0, 2.0);
    const Vector<double> rounding(1.0, 2.0 + 1e-15), error(1.0, 2.0 + 1e-10);
    check(a == b, "exact == on equal vectors");
    check(a != rounding, "exact != on a 1e-15 difference");
    check(ga::approx_equal(a, rounding), "approx_equal absorbs rounding by default");
    check(!ga::approx_equal(a, error), "approx_equal does not absorb a 1e-10 error by default");
    check(ga::approx_equal(a, error, 1e-9), "approx_equal with an explicit, looser tolerance");

    check(Scalar<double>(1.0) == Scalar<double>(1.0), "Scalar ==");
    check(Scalar<double>(1.0) != Scalar<double>(2.0), "Scalar !=");
    check(Bivector<double>(3.0) == Bivector<double>(3.0), "Bivector ==");
    check(Bivector<double>(3.0) != Bivector<double>(-3.0), "Bivector !=");
    check(make_mv(1, 2, 3, 4) == make_mv(1, 2, 3, 4), "Multivector ==");
    check(make_mv(1, 2, 3, 4) != make_mv(1, 2, 3, 5), "Multivector !=");
    check(ga::approx_equal(make_mv(1, 2, 3, 4), make_mv(1, 2, 3, 4 + 1e-15)), "Multivector approx_equal");
    check(ga::approx_equal(Scalar<double>(1.0), Scalar<double>(1.0 + 1e-15)), "Scalar approx_equal");
    check(ga::approx_equal(Bivector<double>(1.0), Bivector<double>(1.0 + 1e-15)), "Bivector approx_equal");

    // A rotor and its negation are the same rotation, but not the same rotor.
    const Rotor<double> r(Scalar<double>(0.6), Bivector<double>(0.8));
    const Rotor<double> neg(Scalar<double>(-0.6), Bivector<double>(-0.8));
    check(r == r, "Rotor ==");
    check(r != neg, "a rotor is not equal to its negation");
    check(!ga::approx_equal(r, neg), "nor approx_equal to it");
    check(ga::approx_equal(r, Rotor<double>(Scalar<double>(0.6 + 1e-15), Bivector<double>(0.8))),
          "Rotor approx_equal");

    // operator<< delegates to to_string() for every type.
    std::ostringstream os;
    os << Scalar<double>(1.5) << '|' << Vector<double>(1, 2) << '|' << Bivector<double>(3)
       << '|' << r << '|' << make_mv(1, 2, 3, 4);
    check(os.str() == "1.5|1*e1 + 2*e2|3*e12|" + r.to_string() + "|(1) + (2*e1 + 3*e2) + (4*e12)",
          "operator<< matches to_string for every type");
}

} // namespace

void run_comparison_and_printing()
{
    test_to_string();
    test_comparison();
}

} // namespace cctest
