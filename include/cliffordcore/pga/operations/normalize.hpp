#pragma once

/**
 * @file cliffordcore/pga/operations/normalize.hpp
 * @brief Unit-size versions of each type.
 *
 * Divides by the Euclidean norm to give a unit version: a plane with a unit
 * normal, a line with a unit direction, a point with weight one, a
 * normalised translator (scalar one). Normalizing a Scalar yields its sign.
 * Ideal elements have no Euclidean size and therefore no normalized form.
 *
 * A Motor is the exception. Its product with its reverse is a study number
 * a + b e0123 rather than a scalar, and dividing by sqrt(a) alone leaves
 * the b part behind, so a motor that has drifted off unit size under
 * repeated composition would still fail m * reverse(m) = 1. normalize(Motor)
 * divides by the square root of the whole study number, which is
 * sqrt(a) (1 + (b / 2a) e0123), and so restores m * reverse(m) = 1 exactly.
 *
 * A zero argument has no direction to preserve, so normalizing one is a
 * precondition violation rather than a value to compute: every overload
 * asserts that the norm is non-zero. The assert compiles out under NDEBUG,
 * so each overload also returns its argument unchanged in that case --
 * a release build must not divide by zero and hand back NaN.
 *
 * There is no normalize(Quadvector): its Euclidean norm is always zero.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <cassert>
#include <cmath>

#include "norm.hpp"
#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"
#include "../translator.hpp"
#include "../motor.hpp"

namespace CliffordCore::PGA
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
     * @brief Normalizes a plane to have a unit normal.
     * @param v The vector to normalize.
     * @return The same plane with a unit normal; its offset e0 scales along.
     * @pre v is not ideal; an ideal plane has no normal to scale and trips an
     *      assert. Under NDEBUG it is returned unchanged instead.
     */
    constexpr Vector<T> normalize(const Vector<T>& v) {
        const T n = norm(v).value;
        assert(n != T(0) && "normalize: the plane is ideal and has no unit normal");
        if (n == T(0)) {
            return v;   // NDEBUG only; the assert above fires otherwise.
        }
        return Vector<T>(v.e0 / n, v.e1 / n, v.e2 / n, v.e3 / n);
    }

    template<typename T>
    /**
     * @brief Normalizes a line to have a unit direction.
     * @param b The bivector to normalize.
     * @return The same line with a unit direction; its moment scales along.
     * @pre b is not ideal; an ideal line has no direction to scale and trips
     *      an assert. Under NDEBUG it is returned unchanged instead.
     */
    constexpr Bivector<T> normalize(const Bivector<T>& b) {
        const T n = norm(b).value;
        assert(n != T(0) && "normalize: the line is ideal and has no unit direction");
        if (n == T(0)) {
            return b;   // NDEBUG only; the assert above fires otherwise.
        }
        return Bivector<T>(b.e01 / n, b.e02 / n, b.e03 / n, b.e12 / n, b.e13 / n, b.e23 / n);
    }

    template<typename T>
    /**
     * @brief Normalizes a point to have weight one.
     * @param t The trivector to normalize.
     * @return The same point with e123 = 1 (or -1 for a negative weight, so
     *         the position components keep their sign convention).
     * @pre t is not ideal; an ideal point has no weight to scale and trips an
     *      assert. Under NDEBUG it is returned unchanged instead.
     */
    constexpr Trivector<T> normalize(const Trivector<T>& t) {
        const T n = norm(t).value;
        assert(n != T(0) && "normalize: the point is ideal and has no unit weight");
        if (n == T(0)) {
            return t;   // NDEBUG only; the assert above fires otherwise.
        }
        return Trivector<T>(t.e012 / n, t.e013 / n, t.e023 / n, t.e123 / n);
    }

    template<typename T>
    /**
     * @brief Normalizes a multivector to unit Euclidean norm.
     * @param m The multivector to normalize.
     * @return The multivector divided by its norm.
     * @pre m has a non-zero Euclidean norm; a purely ideal multivector trips
     *      an assert. Under NDEBUG it is returned unchanged instead.
     */
    constexpr Multivector<T> normalize(const Multivector<T>& m) {
        const T n = norm(m).value;
        assert(n != T(0) && "normalize: the multivector is ideal and has no unit form");
        if (n == T(0)) {
            return m;   // NDEBUG only; the assert above fires otherwise.
        }
        return m / Scalar<T>(n);
    }

    template<typename T>
    /**
     * @brief Normalizes a rotor to unit norm.
     * @param r The rotor to normalize.
     * @return The rotor divided by its norm, so that r * reverse(r) = 1.
     * @pre r is non-zero; a zero rotor trips an assert. Under NDEBUG it is
     *      returned unchanged instead.
     */
    constexpr Rotor<T> normalize(const Rotor<T>& r) {
        const T n = norm(r).value;
        assert(n != T(0) && "normalize: the rotor is zero and has no unit form");
        if (n == T(0)) {
            return r;   // NDEBUG only; the assert above fires otherwise.
        }
        return r / Scalar<T>(n);
    }

    template<typename T>
    /**
     * @brief Normalizes a translator so that its scalar part is one.
     * @param t The translator to normalize.
     * @return The translator divided by its scalar, so that t * reverse(t) = 1.
     * @pre t has a non-zero scalar; a purely ideal translator trips an
     *      assert. Under NDEBUG it is returned unchanged instead.
     */
    constexpr Translator<T> normalize(const Translator<T>& t) {
        const T n = t.scalar.value;
        assert(n != T(0) && "normalize: the translator has a zero scalar and no unit form");
        if (n == T(0)) {
            return t;   // NDEBUG only; the assert above fires otherwise.
        }
        return t / Scalar<T>(n);
    }

    template<typename T>
    /**
     * @brief Normalizes a motor so that m * reverse(m) = 1 exactly.
     * @param m The motor to normalize.
     * @return The motor divided by the square root of its study number
     *         a + b e0123, where a = squared_norm(m) and b is the e0123
     *         part of m * reverse(m). Dividing by sqrt(a) alone is not
     *         enough: it leaves the b term behind.
     * @pre m has a non-zero Euclidean norm; a purely ideal motor trips an
     *      assert. Under NDEBUG it is returned unchanged instead.
     */
    constexpr Motor<T> normalize(const Motor<T>& m) {
        const T a = squared_norm(m).value;
        assert(a != T(0) && "normalize: the motor is ideal and has no unit form");
        if (a == T(0)) {
            return m;   // NDEBUG only; the assert above fires otherwise.
        }
        // The e0123 part of m * reverse(m). Reversing negates the bivector and
        // keeps the scalar and e0123, and the only route to e0123 is the
        // bivector-bivector block of the table plus the scalar-e0123 cross terms.
        const T s = m.scalar.value;
        const T b = T(2) * (s * m.quadvector.e0123
                            - (m.bivector.e01 * m.bivector.e23
                               - m.bivector.e02 * m.bivector.e13
                               + m.bivector.e03 * m.bivector.e12));
        const T inv_sqrt_a = T(1) / std::sqrt(a);
        // 1 / sqrt(a + b I) = (1 / sqrt(a)) (1 - (b / 2a) I), since I^2 = 0.
        const T correction = -b / (T(2) * a);
        const Motor<T> scaled = m * Scalar<T>(inv_sqrt_a);
        // Multiply by (1 + correction * e0123): e0123 commutes with the even
        // subalgebra and only the scalar and Euclidean bivector survive it.
        return Motor<T>(
            scaled.scalar,
            Bivector<T>(
                scaled.bivector.e01 - correction * scaled.bivector.e23,
                scaled.bivector.e02 + correction * scaled.bivector.e13,
                scaled.bivector.e03 - correction * scaled.bivector.e12,
                scaled.bivector.e12,
                scaled.bivector.e13,
                scaled.bivector.e23
            ),
            Quadvector<T>(scaled.quadvector.e0123 + correction * scaled.scalar.value)
        );
    }
} // namespace CliffordCore::PGA
