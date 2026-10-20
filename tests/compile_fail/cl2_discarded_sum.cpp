// expect-error: ignoring return value of
//
// The member operators are [[nodiscard]] too: `a + b;` as a statement computes
// a sum and throws it away, which is never what was meant. See
// cl3_discarded_normalize.cpp for the free-function case and for why the
// pragma below is needed.
//
// check_compile_fail.sh compiles this twice: as written, where it must fail
// with the message above, and with -DCC_T=double, which keeps the result with
// a (void) cast and must succeed.

#include <cliffordcore/cl2.hpp>

#pragma GCC diagnostic error "-Wunused-result"

#ifdef CC_T
#define CC_USE(expr) (void)(expr)
#else
#define CC_USE(expr) expr
#endif

namespace ga = CliffordCore::Cl2;

int main()
{
    const ga::Vector<double> a(1.0, 2.0);
    const ga::Vector<double> b(3.0, 4.0);
    CC_USE(a + b);
}
