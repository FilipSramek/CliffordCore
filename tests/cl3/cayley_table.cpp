// Cl(3,0) tests -- cayley_table: the basis multiplication table, sign by sign.
// The check against an independent blade-product rule runs in properties.cpp,
// in every floating-point type.
//
// One category of the Cl3 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

void test_basis_multiplication_table()
{
    section("basis table");

    // The whole algebra follows from these. Basis: e_i^2 = +1, bivectors stored
    // as (xy, xz, yz) = (e1e2, e1e3, e2e3), pseudoscalar e123 = e1e2e3.
    const Multivector<double> e1  = make_mv(0, 1, 0, 0, 0, 0, 0, 0);
    const Multivector<double> e2  = make_mv(0, 0, 1, 0, 0, 0, 0, 0);
    const Multivector<double> e3  = make_mv(0, 0, 0, 1, 0, 0, 0, 0);
    const Multivector<double> e12 = make_mv(0, 0, 0, 0, 1, 0, 0, 0);
    const Multivector<double> e13 = make_mv(0, 0, 0, 0, 0, 1, 0, 0);
    const Multivector<double> e23 = make_mv(0, 0, 0, 0, 0, 0, 1, 0);
    const Multivector<double> I   = make_mv(0, 0, 0, 0, 0, 0, 0, 1);

    // Vectors square to +1.
    check_close((e1 * e1).scalar.value, 1.0, "e1^2 = +1");
    check_close((e2 * e2).scalar.value, 1.0, "e2^2 = +1");
    check_close((e3 * e3).scalar.value, 1.0, "e3^2 = +1");

    // Bivectors and the pseudoscalar square to -1.
    check_close((e12 * e12).scalar.value, -1.0, "e12^2 = -1");
    check_close((e13 * e13).scalar.value, -1.0, "e13^2 = -1");
    check_close((e23 * e23).scalar.value, -1.0, "e23^2 = -1");
    check_close((I * I).scalar.value, -1.0, "e123^2 = -1");

    // Distinct basis vectors anticommute.
    check_mv(e1 * e2, 0, 0, 0, 0, 1, 0, 0, 0, "e1e2 = e12");
    check_mv(e2 * e1, 0, 0, 0, 0, -1, 0, 0, 0, "e2e1 = -e12");
    check_mv(e1 * e3, 0, 0, 0, 0, 0, 1, 0, 0, "e1e3 = e13");
    check_mv(e2 * e3, 0, 0, 0, 0, 0, 0, 1, 0, "e2e3 = e23");
    check_mv(e1 * e2 * e3, 0, 0, 0, 0, 0, 0, 0, 1, "e1e2e3 = e123");

    // Bivector cross products. These signs are what the (xy, xz, yz) ordering
    // buys, and they are the easiest thing in the library to get backwards.
    check_mv(e12 * e13, 0, 0, 0, 0, 0, 0, -1, 0, "e12 e13 = -e23");
    check_mv(e13 * e12, 0, 0, 0, 0, 0, 0, 1, 0, "e13 e12 = +e23");
    check_mv(e12 * e23, 0, 0, 0, 0, 0, 1, 0, 0, "e12 e23 = +e13");
    check_mv(e23 * e12, 0, 0, 0, 0, 0, -1, 0, 0, "e23 e12 = -e13");
    check_mv(e13 * e23, 0, 0, 0, 0, -1, 0, 0, 0, "e13 e23 = -e12");
    check_mv(e23 * e13, 0, 0, 0, 0, 1, 0, 0, 0, "e23 e13 = +e12");

    // The pseudoscalar is central in 3D: it commutes with everything.
    check_close(mv_difference(e1 * I, I * e1), 0.0, "e123 commutes with e1");
    check_close(mv_difference(e12 * I, I * e12), 0.0, "e123 commutes with e12");
}

} // namespace

void run_cayley_table()
{
    test_basis_multiplication_table();
}

} // namespace cctest
