#pragma once

#include <type_traits>
#include "../scalar.hpp"
#include "../vector3.hpp"
#include "../bivector3.hpp"
#include "../trivector3.hpp"
#include "../multivector3.hpp"
#include "../rotor3.hpp"
#include "dot_product.hpp"
#include "norm.hpp"

namespace CliffordCore
{
    template<typename T>
    /**
     * @brief Computes the inverse of a scalar.
     * @param s The scalar for which to compute the inverse.
     * @return The inverse scalar to the input scalar.
     */
    constexpr Scalar<T> inverse(const Scalar<T>& s) {
        return Scalar<T>(1.0 / s.value);
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a 3D vector.
     * @param v The vector for which to compute the inverse.
     * @return The inverse vector to the input vector.
     */
    constexpr Vector3<T> inverse(const Vector3<T>& v) {
        return v/squared_norm(v);
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a bivector in 3D space.
     * @param b The bivector for which to compute the inverse.
     * @return The inverse bivector to the input bivector.
     */
    constexpr Bivector3<T> inverse(const Bivector3<T>& b) {
        return -b/squared_norm(b);
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a trivector in 3D space.
     * @param t The trivector for which to compute the inverse.
     * @return The inverse trivector to the input trivector.
     */
    constexpr Trivector3<T> inverse(const Trivector3<T>& t) {
        return t/squared_norm(t);
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a multivector in 3D space.
     * @param m The multivector for which to compute the inverse.
     * @return The inverse multivector to the input multivector.
     */
    constexpr Multivector3<T> inverse(const Multivector3<T>& m) {
        return m/squared_norm(m);
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a rotor in 3D space.
     * @param r The rotor for which to compute the inverse.
     * @return The inverse rotor to the input rotor.
     */
    constexpr Rotor3<T> inverse(const Rotor3<T>& r) {
        return r/squared_norm(r);
    }
} // namespace CliffordCore