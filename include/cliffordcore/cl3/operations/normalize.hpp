#pragma once

/**
 * @file cliffordcore/cl3/operations/normalize.hpp
 * @brief Unit-length versions of each type.
 *
 * Divides by the norm to give a unit-length version. Normalizing a Scalar
 * yields its sign. Rotors drift off unit length under repeated composition,
 * faster in float than in double, so renormalise periodically.
 *
 * A zero argument has no direction to preserve, so normalizing one is a
 * precondition violation rather than a value to compute: every overload
 * asserts that the norm is non-zero. The assert compiles out under NDEBUG,
 * so each overload also returns its argument unchanged in that case --
 * a release build must not divide by zero and hand back NaN.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <cassert>

#include "norm.hpp"
#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"

namespace CliffordCore::Cl3
{
    template<typename T>
    /**
     * @brief Normalizes a scalar to unit magnitude.
     * @param s The scalar to normalize.
     * @return +1 or -1, whichever matches the sign of s.
     * @pre s is non-zero; a zero scalar has no sign to preserve and trips an
     *      assert. Under NDEBUG it is returned unchanged instead.
     */
    constexpr Scalar<T> normalize(const Scalar<T>& s) {
        const T n = norm(s).value;
        assert(n != T(0) && "normalize: the scalar is zero and has no unit direction");
        if (n == T(0)) {
            return s;   // NDEBUG only; the assert above fires otherwise.
        }
        return Scalar<T>(s.value / n);
    }

    template<typename T>
    /**
     * @brief Normalizes a 3D vector to have a unit norm.
     * @param v The vector to normalize.
     * @return The normalized vector with a unit norm.
     * @pre v is non-zero; a zero vector has no direction to preserve and trips
     *      an assert. Under NDEBUG it is returned unchanged instead.
     */
    constexpr Vector<T> normalize(const Vector<T>& v) {
        const T n = norm(v).value;
        assert(n != T(0) && "normalize: the vector is zero and has no unit direction");
        if (n == T(0)) {
            return v;   // NDEBUG only; the assert above fires otherwise.
        }
        return Vector<T>(v.x / n, v.y / n, v.z / n);
    }

    template<typename T>
    /**
     * @brief Normalizes a bivector in 3D space to have a unit norm.
     * @param b The bivector to normalize.
     * @return The normalized bivector with a unit norm.
     * @pre b is non-zero; a zero bivector trips an assert. Under NDEBUG it is
     *      returned unchanged instead.
     */
    constexpr Bivector<T> normalize(const Bivector<T>& b) {
        const T n = norm(b).value;
        assert(n != T(0) && "normalize: the bivector is zero and has no unit direction");
        if (n == T(0)) {
            return b;   // NDEBUG only; the assert above fires otherwise.
        }
        return Bivector<T>(b.xy / n, b.xz / n, b.yz / n);
    }

    template<typename T>
    /**
     * @brief Normalizes a trivector in 3D space to have a unit norm.
     * @param t The trivector to normalize.
     * @return The normalized trivector with a unit norm.
     * @pre t is non-zero; a zero trivector trips an assert. Under NDEBUG it is
     *      returned unchanged instead.
     */
    constexpr Trivector<T> normalize(const Trivector<T>& t) {
        const T n = norm(t).value;
        assert(n != T(0) && "normalize: the trivector is zero and has no unit direction");
        if (n == T(0)) {
            return t;   // NDEBUG only; the assert above fires otherwise.
        }
        return Trivector<T>(t.e123 / n);
    }

    template<typename T>
    /**
     * @brief Normalizes a multivector in 3D space to have a unit norm.
     * @param m The multivector to normalize.
     * @return The normalized multivector with a unit norm.
     * @pre m is non-zero; a zero multivector trips an assert. Under NDEBUG it is
     *      returned unchanged instead.
     */
    constexpr Multivector<T> normalize(const Multivector<T>& m) {
        const T n = norm(m).value;
        assert(n != T(0) && "normalize: the multivector is zero and has no unit direction");
        if (n == T(0)) {
            return m;   // NDEBUG only; the assert above fires otherwise.
        }
        return Multivector<T>(
            Scalar<T>(m.scalar.value / n),
            Vector<T>(m.vector.x / n, m.vector.y / n, m.vector.z / n),
            Bivector<T>(m.bivector.xy / n, m.bivector.xz / n, m.bivector.yz / n),
            Trivector<T>(m.trivector.e123 / n)
        );
    }

    template<typename T>
    /**
     * @brief Normalizes a rotor in 3D space to have a unit norm.
     * @param r The rotor to normalize.
     * @return The normalized rotor with a unit norm. Rotors must be unit length
     *         to represent a rotation, so this is the usual way to correct
     *         drift after repeated composition.
     * @pre r is non-zero; a zero rotor trips an assert. Under NDEBUG it is
     *      returned unchanged instead.
     */
    constexpr Rotor<T> normalize(const Rotor<T>& r) {
        const T n = norm(r).value;
        assert(n != T(0) && "normalize: the rotor is zero and has no unit direction");
        if (n == T(0)) {
            return r;   // NDEBUG only; the assert above fires otherwise.
        }
        return Rotor<T>(
            Scalar<T>(r.scalar.value / n),
            Bivector<T>(r.bivector.xy / n, r.bivector.xz / n, r.bivector.yz / n)
        );
    }
} // namespace CliffordCore::Cl3
