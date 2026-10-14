// Cl(2,0) tests -- cayley_table: the 4x4 basis multiplication table: spot
// checks, then all 16 products against an independent implementation of the
// blade-product rule.
//
// One category of the Cl2 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

// ---------------------------------------------------------------------------
// The multiplication table
// ---------------------------------------------------------------------------

void test_basis_multiplication_table()
{
    section("basis multiplication table");

    const Multivector<double> one = basis_blade(0);
    const Multivector<double> e1 = basis_blade(1);
    const Multivector<double> e2 = basis_blade(2);
    const Multivector<double> e12 = basis_blade(3);

    check_mv_close(e1 * e1, one, "e1 e1 == 1");
    check_mv_close(e2 * e2, one, "e2 e2 == 1");
    check_mv_close(e12 * e12, -one, "e12 e12 == -1");
    check_mv_close(e1 * e2, e12, "e1 e2 == e12");
    check_mv_close(e2 * e1, -e12, "e2 e1 == -e12");

    // The pseudoscalar anticommutes with vectors. In Cl(3,0), e123 commutes
    // with everything; this is the line that makes Cl(2,0) its own algebra.
    check_mv_close(e1 * e12, e2, "e1 e12 == +e2");
    check_mv_close(e12 * e1, -e2, "e12 e1 == -e2");
    check_mv_close(e2 * e12, -e1, "e2 e12 == -e1");
    check_mv_close(e12 * e2, e1, "e12 e2 == +e1");
}

void test_cayley_table_against_permutation_rule()
{
    section("cayley table vs permutation rule");

    // All 16 basis products, each checked against the reference rule rather
    // than against a remembered table. This pins every sign in
    // geometric_product(), so a hand-edited term fails the build.
    int nonzero = 0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            int k = -1;
            const int sign = reference_blade_product(i, j, k);
            const Multivector<double> product = basis_blade(i) * basis_blade(j);
            const std::string what = std::string(kBladeNames[i]) + " * " + kBladeNames[j];
            check(k >= 0, what + " lands on a basis blade");
            if (k < 0) continue;
            ++nonzero;
            check_mv_close(product, static_cast<double>(sign) * basis_blade(k),
                           what + " == " + (sign > 0 ? "+" : "-") + kBladeNames[k]);
        }
    }
    // The metric is non-degenerate, so unlike Cl(3,0,1) no product vanishes.
    check(nonzero == 16, "all 16 basis products are non-zero");
}

} // namespace

void run_cayley_table()
{
    test_basis_multiplication_table();
    test_cayley_table_against_permutation_rule();
}

} // namespace cctest
