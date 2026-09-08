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
    constexpr Scalar<T> norm(const Scalar<T>& s) {
        if (s.value == 0) {
            return Scalar<T>(0);
        }

        else {
            return Scalar<T>(1);
        }
    }

    /**
     * @brief Normalizes a 3D vector to have a unit norm.
     * @param v The vector to normalize.
     * @return The normalized vector with a unit norm.
     */
    template<typename T>
    constexpr Vector3<T> normalize(const Vector3<T>& v) {
        auto n = norm(v);
        return Vector3<T>(v.x / n, v.y / n, v.z / n);
    }

    /**
     * @brief Normalizes a bivector in 3D space to have a unit norm.
     * @param b The bivector to normalize.
     * @return The normalized bivector with a unit norm.
     */
    template<typename T>
    constexpr Bivector3<T> normalize(const Bivector3<T>& b) {
        auto n = norm(b);
        return Bivector3<T>(b.xy / n, b.xz / n, b.yz / n);
    }

    /**
     * @brief Normalizes a trivector in 3D space to have a unit norm.
     * @param t The trivector to normalize.
     * @return The normalized trivector with a unit norm.
     */
    template<typename T>
    constexpr Trivector3<T> normalize(const Trivector3<T>& t) {
        auto n = norm(t);
        return Trivector3<T>(t.value / n);
    }

    /**
     * @brief Normalizes a multivector in 3D space to have a unit norm.
     * @param m The multivector to normalize.
     * @return The normalized multivector with a unit norm.
     */
    template<typename T>
    constexpr Multivector3<T> normalize(const Multivector3<T>& m) {
        auto n = norm(m);
        return Multivector3<T>(m.scalar / n, Vector3<T>(m.vector.x / n, m.vector.y / n, m.vector.z / n), Bivector3<T>(m.bivector.xy / n, m.bivector.xz / n, m.bivector.yz / n), Trivector3<T>(m.trivector.value / n));
    }

    /**
     * @brief Normalizes a rotor in 3D space to have a unit norm.
     * @param r The rotor to normalize.
     * @return The normalized rotor with a unit norm.
     */
    template<typename T>
    constexpr Rotor3<T> normalize(const Rotor3<T>& r) {
        auto n = norm(r);
        return Rotor3<T>(r.scalar / n, Bivector3<T>(r.bivector.xy / n, r.bivector.xz / n, r.bivector.yz / n));
    }
} // namespace CliffordCore