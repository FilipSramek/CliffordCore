// expect-error: Vector can only be instantiated with floating-point types
//
// The plainest case: an integer component type is rejected at the type, with a
// message that says which types are accepted.
//
// Compiled twice by check_compile_fail.sh -- see cl2_norm_of_int_vector.cpp.

#include <cliffordcore/cl3.hpp>

#ifndef CC_T
#define CC_T int
#endif

namespace ga = CliffordCore::Cl3;

int main()
{
    const ga::Vector<CC_T> v(1, 2, 3);
    (void)v;
}
