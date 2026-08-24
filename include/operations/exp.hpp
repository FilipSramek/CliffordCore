#pragma once

#include "../rotor3.hpp"
#include "../bivector3.hpp"
#include "../scalar.hpp"
#include "../operations/norm.hpp"

namespace CliffordCore
{
    template<typename T>
    
    constexpr Rotor3<T> exp(const Bivector3<T>& b) {
        Scalar<T> norm_bivector = norm(b);

        if (norm_bivector.value == 0) { // Handle the case when the bivector is zero
            return Rotor3<T>(Scalar<T>(1), Bivector3<T>(0, 0, 0));
        }

        Scalar<T> cos_norm = Scalar<T>(std::cos(norm_bivector.value));

        Bivector3<T> sin_norm = Bivector3<T>(
            std::sin(norm_bivector.value) * (b.xy / norm_bivector.value),
            std::sin(norm_bivector.value) * (b.xz / norm_bivector.value),
            std::sin(norm_bivector.value) * (b.yz / norm_bivector.value)
        );

        return Rotor3<T>(cos_norm, sin_norm);
    }
} // namespace CliffordCore
