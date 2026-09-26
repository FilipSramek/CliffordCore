// expect-error: Scalar can only be instantiated with floating-point types
//
// char is an arithmetic type too, and whether it is signed depends on the
// platform -- one more reason integer-like types are rejected outright.
//
// Compiled twice by check_compile_fail.sh -- see cl2_norm_of_int_vector.cpp.

#include <cliffordcore/cl3.hpp>

#ifndef CC_T
#define CC_T char
#endif

namespace ga = CliffordCore::Cl3;

int main()
{
    const ga::Scalar<CC_T> s(1);
    (void)s;
}
