// Cl(3,0,1) tests -- cayley_table: the 16x16 basis multiplication table: the
// metric and the spot checks from docs/pga.md, then all 256 products against
// an independent implementation of the blade-product rule.
//
// One category of the PGA suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

void test_basis_multiplication_table()
{
    section("basis table");

    const Multivector<double> e0 = basis_blade(1), e1 = basis_blade(2), e2 = basis_blade(3), e3 = basis_blade(4);
    const Multivector<double> zero;

    // The metric: e0 is null, the rest square to +1.
    check(mv_difference(e0 * e0, zero) == 0.0, "e0 * e0 == 0 exactly");
    check_mv_close(e1 * e1, basis_blade(0), "e1 * e1 == 1");
    check_mv_close(e2 * e2, basis_blade(0), "e2 * e2 == 1");
    check_mv_close(e3 * e3, basis_blade(0), "e3 * e3 == 1");

    // Anticommutation of distinct basis vectors, e0 included.
    check_mv_close(e0 * e1, -(e1 * e0), "e0 and e1 anticommute");
    check_mv_close(e1 * e2, -(e2 * e1), "e1 and e2 anticommute");

    // The spot checks from docs/pga.md.
    check_mv_close(basis_blade(5) * basis_blade(8), basis_blade(6), "e01 * e12 == +e02");
    check_mv_close(basis_blade(8) * basis_blade(5), -basis_blade(6), "e12 * e01 == -e02");
    check_mv_close(basis_blade(14) * basis_blade(15), basis_blade(1), "e123 * e0123 == +e0");
    check_mv_close(basis_blade(15) * basis_blade(14), -basis_blade(1), "e0123 * e123 == -e0");
    check_mv_close(basis_blade(11) * basis_blade(14), -basis_blade(7), "e012 * e123 == -e03");
    check(mv_difference(basis_blade(11) * basis_blade(12), zero) == 0.0, "e012 * e013 == 0 exactly");
    check(mv_difference(basis_blade(15) * basis_blade(15), zero) == 0.0, "e0123 * e0123 == 0 exactly");
    check_mv_close(basis_blade(14) * basis_blade(14), -basis_blade(0), "e123 * e123 == -1");
    check_mv_close(basis_blade(8) * basis_blade(15), -basis_blade(7), "e12 * e0123 == -e03");
    check_mv_close(basis_blade(9) * basis_blade(15), basis_blade(6), "e13 * e0123 == +e02");
    check_mv_close(basis_blade(10) * basis_blade(15), -basis_blade(5), "e23 * e0123 == -e01");

    // The pseudoscalar commutes with every bivector and anticommutes with every vector.
    const Multivector<double> I = basis_blade(15);
    for (int b = 5; b <= 10; ++b) {
        check_mv_close(I * basis_blade(b), basis_blade(b) * I, std::string("e0123 commutes with ") + kBladeNames[b]);
    }
    for (int v = 1; v <= 4; ++v) {
        check_mv_close(I * basis_blade(v), -(basis_blade(v) * I), std::string("e0123 anticommutes with ") + kBladeNames[v]);
    }
}

void test_cayley_table_against_permutation_rule()
{
    section("Cayley table vs permutation rule");

    // Every one of the 256 basis products, against the independent reference.
    // One check per entry keeps a failure pointing at the exact term.
    int nonzero = 0;
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            int k = -1;
            const int sign = reference_blade_product(i, j, k);
            const Multivector<double> product = basis_blade(i) * basis_blade(j);
            const std::string what = std::string(kBladeNames[i]) + " * " + kBladeNames[j];
            if (sign == 0) {
                check(mv_difference(product, Multivector<double>()) == 0.0, what + " == 0");
            } else {
                ++nonzero;
                check_mv_close(product, static_cast<double>(sign) * basis_blade(k),
                               what + " == " + (sign > 0 ? "+" : "-") + kBladeNames[k]);
            }
        }
    }
    check(nonzero == 192, "192 of the 256 basis products are non-zero");
}

} // namespace

void run_cayley_table()
{
    test_basis_multiplication_table();
    test_cayley_table_against_permutation_rule();
}

} // namespace cctest
