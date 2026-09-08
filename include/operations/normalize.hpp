#pragma once

#include "norm.hpp"
#include "../scalar.hpp"
#include "../vector3.hpp"
#include "../bivector3.hpp"
#include "../trivector3.hpp"
#include "../multivector3.hpp"
#include "../rotor3.hpp"

namespace CliffordCore
{
    template<typename T>
    /**
     * @brief Normalizes a 3D vector to have a unit norm.
     * @param v The vector to normalize.
     * @return The normalized vector with a unit norm. A zero vector is returned
     *         unchanged, since it has no direction to preserve.
     */
    constexpr Vector3<T> normalize(const Vector3<T>& v) {
        const T n = norm(v).value;
        if (n == T(0)) {
            return v;
        }
        return Vector3<T>(v.x / n, v.y / n, v.z / n);
    }

    template<typename T>
    /**
     * @brief Normalizes a bivector in 3D space to have a unit norm.
     * @param b The bivector to normalize.
     * @return The normalized bivector with a unit norm. A zero bivector is
     *         returned unchanged.
     */
    constexpr Bivector3<T> normalize(const Bivector3<T>& b) {
        const T n = norm(b).value;
        if (n == T(0)) {
            return b;
        }
        return Bivector3<T>(b.xy / n, b.xz / n, b.yz / n);
    }

    template<typename T>
    /**
     * @brief Normalizes a trivector in 3D space to have a unit norm.
     * @param t The trivector to normalize.
     * @return The normalized trivector with a unit norm. A zero trivector is
     *         returned unchanged.
     */
    constexpr Trivector3<T> normalize(const Trivector3<T>& t) {
        const T n = norm(t).value;
        if (n == T(0)) {
            return t;
        }
        return Trivector3<T>(t.e123 / n);
    }

    template<typename T>
    /**
     * @brief Normalizes a multivector in 3D space to have a unit norm.
     * @param m The multivector to normalize.
     * @return The normalized multivector with a unit norm. A zero multivector is
     *         returned unchanged.
     */
    constexpr Multivector3<T> normalize(const Multivector3<T>& m) {
        const T n = norm(m).value;
        if (n == T(0)) {
            return m;
        }
        return Multivector3<T>(
            Scalar<T>(m.scalar.value / n),
            Vector3<T>(m.vector.x / n, m.vector.y / n, m.vector.z / n),
            Bivector3<T>(m.bivector.xy / n, m.bivector.xz / n, m.bivector.yz / n),
            Trivector3<T>(m.trivector.e123 / n)
        );
    }

    template<typename T>
    /**
     * @brief Normalizes a rotor in 3D space to have a unit norm.
     * @param r The rotor to normalize.
     * @return The normalized rotor with a unit norm. A zero rotor is returned
     *         unchanged. Rotors must be unit length to represent a rotation, so
     *         this is the usual way to correct drift after repeated composition.
     */
    constexpr Rotor3<T> normalize(const Rotor3<T>& r) {
        const T n = norm(r).value;
        if (n == T(0)) {
            return r;
        }
        return Rotor3<T>(
            Scalar<T>(r.scalar.value / n),
            Bivector3<T>(r.bivector.xy / n, r.bivector.xz / n, r.bivector.yz / n)
        );
    }
} // namespace CliffordCore
