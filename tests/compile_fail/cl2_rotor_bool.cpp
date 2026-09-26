// expect-error: Rotor can only be instantiated with floating-point types
//
// bool is an arithmetic type, so the old std::is_arithmetic guard let it
// through, and the whole API instantiated for it: sums of two `true`
// components were stored straight back into a bool.
//
// Compiled twice by check_compile_fail.sh -- see cl2_norm_of_int_vector.cpp.

#include <cliffordcore/cl2.hpp>

#ifndef CC_T
#define CC_T bool
#endif

namespace ga = CliffordCore::Cl2;

int main()
{
    const ga::Rotor<CC_T> r;
    (void)r;
}
