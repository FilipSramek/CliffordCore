#pragma once

#include "../rotor3.hpp"
#include "../bivector3.hpp"
#include "../scalar.hpp"
#include "../operations/norm.hpp"

namespace CliffordCore
{
    template<typename T>
    
    constexpr Bivector3<T> log(const Rotor3<T>& r) {
        Scalar<T> arccos_scalar = Scalar<T> (std::acos(r.scalar.value));

        Scalar<T> Magnitude_bivector = Bivector3<T>(r.bivector.xy, r.bivector.xz, r.bivector.yz).magnitude();

        if (Magnitude_bivector.value == 0) { // Handle the case when the bivector is zero
            return Bivector3<T>(0, 0, 0);
        }

        log = Bivector3<T>(
            (arccos_scalar.value / Magnitude_bivector.value) * r.bivector.xy,
            (arccos_scalar.value / Magnitude_bivector.value) * r.bivector.xz,
            (arccos_scalar.value / Magnitude_bivector.value) * r.bivector.yz
        );

        return log;
    }
} // namespace CliffordCore
