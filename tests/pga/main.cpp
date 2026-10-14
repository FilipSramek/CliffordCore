// CliffordCore test suite for Cl(3,0,1), namespace CliffordCore::PGA.
//
// Build and run (Git Bash):
//   g++ -std=c++17 -Iinclude tests/pga/*.cpp -o build/test_pga && ./build/test_pga.exe
//   ...or just ./build.sh, or ./build.sh pga for this suite alone
//
// Structure: one file per kind of test, each its own translation unit, all
// linked into this one binary. support.hpp includes every PGA header and
// holds the shared typed helpers; the assertion harness is
// tests/common/harness.hpp, shared by every suite but including no algebra.
// Every check records a pass/fail and the run continues, so a single build
// reports every problem rather than stopping at the first.
//
// The sign conventions asserted in this suite are the ones documented in docs/pga.md;
// every one of them was verified against an independent model of the algebra
// before being written down, so a failure means the code, not the expectation.

#include "support.hpp"

int main()
{
    using namespace cctest;

    run_types();
    run_arithmetic();
    run_products();
    run_cayley_table();
    run_incidence();
    run_grades_and_contractions();
    run_unary_ops();
    run_motions();
    run_geometry();
    run_comparison_and_printing();
    run_properties();
    run_instantiation_sweep();

    std::cout << g_checks << " checks, " << g_failures << " failed.\n";
    if (g_failures != 0) {
        std::cout << "TESTS FAILED\n";
        return 1;
    }
    std::cout << "All CliffordCore PGA tests passed\n";
    return 0;
}
