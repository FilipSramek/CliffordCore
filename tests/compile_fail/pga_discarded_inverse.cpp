// expect-error: ignoring return value of
//
// inverse() returns the inverse motor and leaves its argument alone; discarding
// the result is a mistake the compiler can see. See
// cl3_discarded_normalize.cpp for why the pragma below is needed.
//
// check_compile_fail.sh compiles this twice: as written, where it must fail
// with the message above, and with -DCC_T=double, which keeps the result with
// a (void) cast and must succeed.

#include <cliffordcore/pga.hpp>

#pragma GCC diagnostic error "-Wunused-result"

#ifdef CC_T
#define CC_USE(expr) (void)(expr)
#else
#define CC_USE(expr) expr
#endif

namespace ga = CliffordCore::PGA;

int main()
{
    const ga::Motor<double> m = ga::exp(ga::Bivector<double>(0.1, 0.2, 0.3, 0.4, 0.5, 0.6));
    CC_USE(ga::inverse(m));
}
