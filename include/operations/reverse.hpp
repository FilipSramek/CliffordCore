#pragma once

#include "../vector3.hpp"
#include "../bivector3.hpp"

namespace CliffordCore
{
    template<typename T>
    /**
     * @brief Computes the geometric reverse of a 3D vector.
     * @param v The vector for which to compute the reverse.
     * @return The reverse of the vector v. (reverse(v) = v)
     */
    constexpr Vector3<T> reverse(const Vector3<T>& v) {
        return Vector3<T>(v.x, v.y, v.z);
    }

    template<typename T>
    /**
     * @brief Computes the geometric reverse of a bivector in 3D space.
     * @param b The bivector for which to compute the reverse.
     * @return The reverse of the bivector b. (reverse(b) = -b)
     */
    constexpr Bivector3<T> reverse(const Bivector3<T>& b) {
        return Bivector3<T>(-b.xy, -b.xz, -b.yz);
    }

    template<typename T>
    /**
     * @brief Computes the geometric reverse of a rotor in 3D space.
     * @param r The rotor for which to compute the reverse.
     * @return The reverse of the rotor r. (reverse(r) = scalar(r) - bivector(r))
     */

    constexpr Rotor3<T> reverse(const Rotor3<T>& r) {
        return Rotor3<T>(r.scalar, Bivector3<T>(-r.bivector.xy, -r.bivector.xz, -r.bivector.yz));
    }
} // namespace CliffordCore