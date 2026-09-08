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
        const T s  = r.scalar.value;
        const T xy = r.bivector.xy;
        const T xz = r.bivector.xz;
        const T yz = r.bivector.yz;

        // Expanding R v ~R. Note which axis each bivector turns: e12 (xy) spans
        // the xy plane, so it rotates about z -- not about x. Pairing xy with
        // the x axis is the easy mistake here.
        return Vector3<T>(
            (s*s - xy*xy - xz*xz + yz*yz) * v.x
                + 2 * (s*xy - xz*yz) * v.y
                + 2 * (s*xz + xy*yz) * v.z,

            2 * (-s*xy - xz*yz) * v.x
                + (s*s - xy*xy + xz*xz - yz*yz) * v.y
                + 2 * (s*yz - xy*xz) * v.z,

            2 * (-s*xz + xy*yz) * v.x
                + 2 * (-s*yz - xy*xz) * v.y
                + (s*s + xy*xy - xz*xz - yz*yz) * v.z
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