// expect-error: Vector can only be instantiated with floating-point types
//
// The case that motivated the floating-point-only contract. Before it,
// norm(Vector<int>(1, 1)) compiled and returned 1 -- sqrt(2) truncated -- with
// no diagnostic. A free function needs no guard of its own: naming Vector<int>
// instantiates the class, and the class guard fires.
//
// check_compile_fail.sh compiles this twice: as written, where it must fail
// with the message above, and with -DCC_T=double, where it must succeed. The
// second run proves the failure is caused by the type and nothing else.

#include <cliffordcore/cl2.hpp>

#ifndef CC_T
#define CC_T int
#endif

namespace ga = CliffordCore::Cl2;

int main()
{
    const auto n = ga::norm(ga::Vector<CC_T>(1, 1));
    (void)n;
}
