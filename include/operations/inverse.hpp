#pragma once

/**
 * @file inverse.hpp
 * @brief Multiplicative inverses, and Clifford conjugation.
 *
 * Inverses carry the minus signs the squares imply: v^-1 = v/|v|^2, but
 * B^-1 = -B/|B|^2 and t^-1 = -t/|t|^2. A Multivector3 inverts through
 * Clifford conjugation rather than by dividing by its norm.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <type_traits>
#include "../scalar.hpp"
#include "../vector3.hpp"
#include "../bivector3.hpp"
#include "../trivector3.hpp"
#include "../multivector3.hpp"
#include "../rotor3.hpp"
#include "dot_product.hpp"
#include "norm.hpp"
#include "reverse.hpp"
#include "involutions.hpp"
#include "geometric_product.hpp"

namespace CliffordCore
{
    template<typename T>
    /**
     * @brief Computes the inverse of a scalar.
     * @param s The scalar for which to compute the inverse.
     * @return The inverse scalar to the input scalar.
     */
    constexpr Scalar<T> inverse(const Scalar<T>& s) {
        return Scalar<T>(T(1) / s.value);
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
     * @return The inverse trivector to the input trivector. The pseudoscalar
     *         squares to -1, so the inverse carries a minus sign: t^-1 = -t/|t|^2.
     */
    constexpr Trivector3<T> inverse(const Trivector3<T>& t) {
        return -t/squared_norm(t);
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a multivector in 3D space.
     * @param m The multivector for which to compute the inverse.
     * @return The inverse multivector, satisfying m * inverse(m) = 1.
     *
     * Dividing by the squared norm only works for a single grade. In general,
     * m * conjugate(m) collapses to a scalar plus a pseudoscalar, a + b*e123.
     * Because e123 squares to -1 that behaves like a complex number, so it
     * inverts as (a - b*e123)/(a^2 + b^2), and
     * inverse(m) = conjugate(m) * inverse(m * conjugate(m)).
     *
     * The multivector is not invertible when a^2 + b^2 is zero; this returns a
     * zero multivector in that case rather than dividing by zero.
     */
    constexpr Multivector3<T> inverse(const Multivector3<T>& m) {
        const Multivector3<T> conj = conjugate(m);
        const Multivector3<T> collapsed = m * conj;

        const T a = collapsed.scalar.value;
        const T b = collapsed.trivector.e123;
        const T denominator = a * a + b * b;

        if (denominator == T(0)) {
            return Multivector3<T>();
        }

        const Multivector3<T> collapsedInverse(
            Scalar<T>(a / denominator),
            Vector3<T>(),
            Bivector3<T>(),
            Trivector3<T>(-b / denominator)
        );

        return conj * collapsedInverse;
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a rotor in 3D space.
     * @param r The rotor for which to compute the inverse.
     * @return The inverse rotor to the input rotor: r^-1 = reverse(r)/|r|^2,
     *         since r * reverse(r) = |r|^2. For a unit rotor this is just the
     *         reverse, which is why rotations undo by reversing.
     */
    constexpr Rotor3<T> inverse(const Rotor3<T>& r) {
        return reverse(r)/squared_norm(r);
    }
} // namespace CliffordCore