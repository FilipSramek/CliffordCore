// CliffordCore test suite for Cl(2,0), namespace CliffordCore::Cl2.
//
// Build and run (Git Bash):
//   g++ -std=c++17 -Iinclude tests/cl2/*.cpp -o build/test_cl2 && ./build/test_cl2.exe
//   ...or just ./build.sh, or ./build.sh cl2 for this suite alone
//
// Structure: one file per kind of test, each its own translation unit, all
// linked into this one binary. support.hpp includes every Cl2 header and
// holds the shared typed helpers; the assertion harness is
// tests/common/harness.hpp, shared by every suite but including no algebra.
// Every check records a pass/fail and the run continues, so a single build
// reports every problem rather than stopping at the first.
//
// The sign conventions asserted in this suite are the ones documented in docs/cl2.md;
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
    run_grades_and_contractions();
    run_unary_ops();
    run_rotations();
    run_geometry();
    run_complex_isomorphism();
    run_comparison_and_printing();
    run_properties();
    run_special_values();
    run_instantiation_sweep();

    std::cout << g_checks << " checks, " << g_failures << " failed.\n";
    if (g_failures != 0) {
        std::cout << "TESTS FAILED\n";
        return 1;
    }
    std::cout << "All CliffordCore Cl2 tests passed\n";
    return 0;
}
