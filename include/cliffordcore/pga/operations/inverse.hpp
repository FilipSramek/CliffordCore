#pragma once

/**
 * @file cliffordcore/pga/operations/inverse.hpp
 * @brief Multiplicative inverses in the degenerate metric.
 *
 * Not everything in Cl(3,0,1) has an inverse: e0 and everything ideal
 * squares to zero, and the pseudoscalar e0123 with it. The invertible
 * things invert like this:
 *
 *   - Vector, Trivector, Rotor, Translator: x ~x is a scalar, so
 *     x^-1 = reverse(x) / squared_norm(x), exactly as in Cl3.
 *   - Bivector and Motor: x ~x is a study number a + b e0123, and since
 *     e0123^2 = 0 its inverse is (a - b e0123) / a^2, so
 *     x^-1 = reverse(x) (a - b e0123) / a^2. Dropping the b term gives a
 *     wrong answer for any motor that is not already normalised.
 *   - Multivector: the general n = 4 formula of Hitzer and Sangwine, which
 *     survives the degenerate metric. With N = m * conjugate(m) (grades 0,
 *     3 and 4 only) and N' = N with grades 3 and 4 negated, N N' is a pure
 *     scalar and m^-1 = conjugate(m) N' / <N N'>_0.
 *
 * A non-invertible argument returns zero rather than dividing by zero, as
 * Cl3's multivector inverse does. There is no inverse(Quadvector).
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <type_traits>
#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../quadvector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"
#include "../translator.hpp"
#include "../motor.hpp"
#include "norm.hpp"
#include "reverse.hpp"
#include "involutions.hpp"
#include "geometric_product.hpp"

namespace CliffordCore::PGA
{
    template<typename T>
    /**
     * @brief Computes the inverse of a scalar.
     * @param s The scalar for which to compute the inverse.
     * @return 1 / s.
     */
    constexpr Scalar<T> inverse(const Scalar<T>& s) {
        return Scalar<T>(T(1) / s.value);
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a plane.
     * @param v The vector for which to compute the inverse.
     * @return v / squared_norm(v), since a plane squares to the squared
     *         length of its normal. An ideal plane has no inverse; this
     *         returns zero for one.
     */
    constexpr Vector<T> inverse(const Vector<T>& v) {
        const T a = squared_norm(v).value;
        if (a == T(0)) {
            return Vector<T>();
        }
        return v / Scalar<T>(a);
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a line.
     * @param b The bivector for which to compute the inverse.
     * @return reverse(b) (a - c e0123) / a^2 where a + c e0123 = b * reverse(b).
     *         For a line (a simple bivector) c is zero and this is -b / |b|^2,
     *         as in Cl3; for a general bivector the e0123 correction matters.
     *         An ideal line has no inverse; this returns zero for one.
     */
    constexpr Bivector<T> inverse(const Bivector<T>& b) {
        const T a = squared_norm(b).value;
        if (a == T(0)) {
            return Bivector<T>();
        }
        // b * reverse(b) = -b^2 = a - 2 (e01 e23 - e02 e13 + e03 e12) e0123.
        const T c = -T(2) * (b.e01 * b.e23 - b.e02 * b.e13 + b.e03 * b.e12);
        // reverse(b) (a - c e0123): e0123 turns the Euclidean part into an
        // ideal one via e12 e0123 = -e03, e13 e0123 = +e02, e23 e0123 = -e01.
        return Bivector<T>(
            (-a * b.e01 - c * b.e23) / (a * a),
            (-a * b.e02 + c * b.e13) / (a * a),
            (-a * b.e03 - c * b.e12) / (a * a),
            -b.e12 / a,
            -b.e13 / a,
            -b.e23 / a
        );
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a point.
     * @param t The trivector for which to compute the inverse.
     * @return reverse(t) / squared_norm(t) = -t / e123^2, since a point
     *         times its reverse is its squared weight. An ideal point has no
     *         inverse; this returns zero for one.
     */
    constexpr Trivector<T> inverse(const Trivector<T>& t) {
        const T a = squared_norm(t).value;
        if (a == T(0)) {
            return Trivector<T>();
        }
        return reverse(t) / Scalar<T>(a);
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a multivector.
     * @param m The multivector for which to compute the inverse.
     * @return The inverse multivector, satisfying m * inverse(m) = 1, or zero
     *         if m is not invertible.
     *
     * Hitzer and Sangwine's formula for four dimensions: N = m * conjugate(m)
     * has only grades 0, 3 and 4; negating grades 3 and 4 gives N', and N N'
     * is a pure scalar. Then inverse(m) = conjugate(m) N' / <N N'>_0. The
     * derivation uses only the grading and the involutions, not the metric,
     * so it holds with e0^2 = 0 -- the tests check it on full 16-component
     * values. The denominator is zero exactly when m has no inverse (e0, for
     * instance), and zero is returned then rather than NaN.
     */
    constexpr Multivector<T> inverse(const Multivector<T>& m) {
        const Multivector<T> conj = conjugate(m);
        const Multivector<T> n = m * conj;
        const Multivector<T> nPrime(n.scalar, n.vector, n.bivector, -n.trivector, -n.quadvector);
        const T denominator = (n * nPrime).scalar.value;

        if (denominator == T(0)) {
            return Multivector<T>();
        }

        return (conj * nPrime) / Scalar<T>(denominator);
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a rotor.
     * @param r The rotor for which to compute the inverse.
     * @return reverse(r) / squared_norm(r), since r * reverse(r) = |r|^2. For
     *         a unit rotor this is just the reverse.
     */
    constexpr Rotor<T> inverse(const Rotor<T>& r) {
        const T a = squared_norm(r).value;
        if (a == T(0)) {
            return Rotor<T>();
        }
        return reverse(r) / Scalar<T>(a);
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a translator.
     * @param t The translator for which to compute the inverse.
     * @return reverse(t) / scalar^2, since t * reverse(t) = scalar^2 (the
     *         ideal part squares to zero). For a normalised translator this
     *         is just the reverse: the opposite translation. A purely ideal
     *         translator has no inverse; this returns zero for one.
     */
    constexpr Translator<T> inverse(const Translator<T>& t) {
        const T a = squared_norm(t).value;
        if (a == T(0)) {
            return Translator<T>();
        }
        return reverse(t) / Scalar<T>(a);
    }

    template<typename T>
    /**
     * @brief Computes the inverse of a motor.
     * @param m The motor for which to compute the inverse.
     * @return reverse(m) (a - c e0123) / a^2 where a + c e0123 = m * reverse(m)
     *         is the motor's study number. For a unit motor this is just the
     *         reverse; for one that has drifted, the e0123 correction is what
     *         keeps m * inverse(m) = 1. A purely ideal motor has no inverse;
     *         this returns zero for one.
     */
    constexpr Motor<T> inverse(const Motor<T>& m) {
        const T a = squared_norm(m).value;
        if (a == T(0)) {
            return Motor<T>();
        }
        const T s = m.scalar.value;
        const T p = m.quadvector.e0123;
        // m * reverse(m) = (s^2 + |B_e|^2) + 2 (s p - (e01 e23 - e02 e13 + e03 e12)) e0123.
        const T c = T(2) * (s * p - (m.bivector.e01 * m.bivector.e23
                                     - m.bivector.e02 * m.bivector.e13
                                     + m.bivector.e03 * m.bivector.e12));
        // reverse(m) (a - c e0123) = a s - a B + c (B e0123) + (a p - c s) e0123,
        // with B e0123 = B_e e0123 mapping (e12, e13, e23) to (-e03, +e02, -e01).
        return Motor<T>(
            Scalar<T>(s / a),
            Bivector<T>(
                (-a * m.bivector.e01 - c * m.bivector.e23) / (a * a),
                (-a * m.bivector.e02 + c * m.bivector.e13) / (a * a),
                (-a * m.bivector.e03 - c * m.bivector.e12) / (a * a),
                -m.bivector.e12 / a,
                -m.bivector.e13 / a,
                -m.bivector.e23 / a
            ),
            Quadvector<T>((a * p - c * s) / (a * a))
        );
    }
} // namespace CliffordCore::PGA
