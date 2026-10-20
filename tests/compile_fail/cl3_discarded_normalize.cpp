// expect-error: ignoring return value of
//
// normalize() returns the unit rotor; it does not change its argument. Called
// as a statement, by someone expecting it to normalize in place, it used to
// compile silently and leave the rotor as it was. Every function that returns
// a value is [[nodiscard]], which turns that into a warning, and the pragma
// below makes the warning an error so this check can see it.
//
// check_compile_fail.sh compiles this twice: as written, where it must fail
// with the message above, and with -DCC_T=double, which keeps the result with
// a (void) cast and must succeed. The second run proves the failure is caused
// by the discarded result and nothing else.

#include <cliffordcore/cl3.hpp>

#pragma GCC diagnostic error "-Wunused-result"

#ifdef CC_T
#define CC_USE(expr) (void)(expr)
#else
#define CC_USE(expr) expr
#endif

namespace ga = CliffordCore::Cl3;

int main()
{
    const ga::Rotor<double> r(ga::Scalar<double>(2.0), ga::Bivector<double>(0.0, 0.0, 0.0));
    CC_USE(ga::normalize(r));
}
