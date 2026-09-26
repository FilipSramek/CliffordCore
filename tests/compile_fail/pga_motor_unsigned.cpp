// expect-error: Motor can only be instantiated with floating-point types
//
// Under the old std::is_arithmetic guard, unsigned types got further than any
// other integer type before failing -- on an ambiguous std::abs deep inside
// norm.hpp. Now they stop at the type, with a message that says why.
//
// Compiled twice by check_compile_fail.sh -- see cl2_norm_of_int_vector.cpp.

#include <cliffordcore/pga.hpp>

#ifndef CC_T
#define CC_T unsigned
#endif

namespace ga = CliffordCore::PGA;

int main()
{
    const ga::Motor<CC_T> m;
    (void)m;
}
