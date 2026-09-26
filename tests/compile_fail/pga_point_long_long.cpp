// expect-error: Trivector can only be instantiated with floating-point types
//
// The geometric aliases carry the guard with them: Point<T> is Trivector<T>,
// so the diagnostic names Trivector.
//
// Compiled twice by check_compile_fail.sh -- see cl2_norm_of_int_vector.cpp.

#include <cliffordcore/pga.hpp>

#ifndef CC_T
#define CC_T long long
#endif

namespace ga = CliffordCore::PGA;

int main()
{
    const ga::Point<CC_T> p;
    (void)p;
}
