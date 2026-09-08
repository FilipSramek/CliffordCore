#pragma once

#include <type_traits>
#include "../vector3.hpp"
#include "../multivector3.hpp"
#include "../rotor3.hpp"
#include "../bivector3.hpp"
#include "../scalar.hpp" 
#include "dot_product.hpp"
#include "wedge_product.hpp"

namespace CliffordCore
{
    template<typename T>
    /**
     * @brief Computes the geometric product of two 3D vectors, resulting in a multivector.
     * @param a The first vector.
     * @param b The second vector.
     * @return The resulting multivector from the geometric product = dot_product(a, b) + wedge_product(a, b) so that the scalar part is the dot product and the bivector part is the wedge product.
     */
    constexpr Multivector3<T> geometric_product(const Vector3<T>& a, const Vector3<T>& b) {
        return Multivector3<T>(
            a | b,
            Vector3<T>(0, 0, 0),
            a ^ b,
            Trivector3<T>(0)
        );
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two 3D vectors, resulting in a rotor.
     * @param a The first vector.
     * @param b The second vector.
     * @return The resulting rotor from the geometric product.
     */
    constexpr Rotor3<T> operator*(const Vector3<T>& a, const Vector3<T>& b) {
        return Rotor3<T>(
            a | b,
            a ^ b
        );
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two 3D vectors, resulting in a multivector.
     * @param a The first vector.
     * @param b The second vector.
     * @return The resulting multivector from the geometric product = dot_product(a, b) + wedge_product(a, b) so that the scalar part is the dot product and the bivector part is the wedge product.
     */
    constexpr Multivector3<T> operator*(const Vector3<T>& a, const Vector3<T>& b) {
        return Multivector3<T>(
            a | b,
            Vector3<T>(0, 0, 0),
            a ^ b,
            Trivector3<T>(0)
        );
    }

    template<typename T>
    /**
     * @brief Computes the geometric product of two 3D rotors, resulting in a 3D rotor. Analogous to two conecutive rotations.
     * @param r The first rotor.
     * @param s The second rotor.
     * @return The resulting rotor from the geometric product.
     */
    constexpr Rotor3<T> operator*(const Rotor3<T>& r, const Rotor3<T>& s) {

    Scalar<T> scalar =
        r.scalar * s.scalar
        - (r.bivector.xy * s.bivector.xy
         + r.bivector.xz * s.bivector.xz
         + r.bivector.yz * s.bivector.yz);

    Bivector3<T> bivector(
        r.scalar * s.bivector.xy
        + r.bivector.xy * s.scalar
        - r.bivector.xz * s.bivector.yz
        + r.bivector.yz * s.bivector.xz,

        r.scalar * s.bivector.xz
        + r.bivector.xz * s.scalar
        + r.bivector.xy * s.bivector.yz
        - r.bivector.yz * s.bivector.xy,

        r.scalar * s.bivector.yz
        + r.bivector.yz * s.scalar
        - r.bivector.xy * s.bivector.xz
        + r.bivector.xz * s.bivector.xy
    );

        return Rotor3<T>(scalar, bivector);
    }
} // namespace CliffordCore