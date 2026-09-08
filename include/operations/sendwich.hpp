#pragma once

#include <type_traits>
#include "../vector3.hpp"
#include "../rotor3.hpp"

namespace CliffordCore
{
    template<typename T>
    /**
     * @brief Computes the sandwich product of a 3D vector and a 3D rotor (R*v*Reverse(R)).
     * @param v The vector to be transformed.
     * @param r The rotor to perform the transformation with.
     * @return The resulting vector after the sandwich product.
     */
    constexpr Vector3<T> sandwich(const Vector3<T>& v, const Rotor3<T>& r) {
        const T s  = r.scalar;
        const T xy = r.bivector.xy;
        const T xz = r.bivector.xz;
        const T yz = r.bivector.yz;

        return Vector3<T>(
            (s*s + xy*xy - xz*xz - yz*yz) * v.x
                + 2 * (xy*xz - s*yz) * v.y
                + 2 * (xy*yz + s*xz) * v.z,

            2 * (xy*xz + s*yz) * v.x
                + (s*s - xy*xy + xz*xz - yz*yz) * v.y
                + 2 * (xz*yz - s*xy) * v.z,

            2 * (xy*yz - s*xz) * v.x
                + 2 * (xz*yz + s*xy) * v.y
                + (s*s - xy*xy - xz*xz + yz*yz) * v.z
        );
    }

    /**
     * @brief Rotates a 3D vector using a 3D rotor.
     * @param v The vector to be rotated.
     * @param r The rotor to perform the rotation with.
     * @return The resulting vector after rotation.
     */
    template<typename T>
    constexpr Vector3<T> rotate(const Vector3<T>& v, const Rotor3<T>& r) {
        return sandwich(v, r);
    }
} // namespace CliffordCore