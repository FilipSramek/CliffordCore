#pragma once

/**
 * @file involutions.hpp
 * @brief The grade involution and Clifford conjugation.
 */

#include "../scalar.hpp"
#include "../vector3.hpp"
#include "../bivector3.hpp"
#include "../trivector3.hpp"
#include "../multivector3.hpp"
#include "../rotor3.hpp"
#include "reverse.hpp"

// Cl(3,0) has three sign-flipping involutions, and they come as a set. Each
// scales grade k by a sign that depends only on k:
//
//              grade:            0    1    2    3
//   reverse    (-1)^(k(k-1)/2)   +    +    -    -
//   involute   (-1)^k            +    -    +    -
//   conjugate  (-1)^(k(k+1)/2)   +    -    -    +
//
// So conjugate == reverse composed with involute, in either order. The test
// suite asserts that identity rather than trusting the table.
//
// reverse() lives in reverse.hpp, which this header includes, so including
// this one gives you all three.

namespace CliffordCore
{
    template<typename T>
    /**
     * @brief Applies the grade involution to a scalar.
     * @param s The scalar to involve.
     * @return The scalar unchanged; grade 0 keeps its sign.
     */
    constexpr Scalar<T> involute(const Scalar<T>& s) {
        return s;
    }

    template<typename T>
    /**
     * @brief Applies the grade involution to a vector.
     * @param v The vector to involve.
     * @return The negated vector; grade 1 flips sign.
     */
    constexpr Vector3<T> involute(const Vector3<T>& v) {
        return Vector3<T>(-v.x, -v.y, -v.z);
    }

    template<typename T>
    /**
     * @brief Applies the grade involution to a bivector.
     * @param b The bivector to involve.
     * @return The bivector unchanged; grade 2 keeps its sign.
     */
    constexpr Bivector3<T> involute(const Bivector3<T>& b) {
        return b;
    }

    template<typename T>
    /**
     * @brief Applies the grade involution to a trivector.
     * @param t The trivector to involve.
     * @return The negated trivector; grade 3 flips sign.
     */
    constexpr Trivector3<T> involute(const Trivector3<T>& t) {
        return Trivector3<T>(-t.e123);
    }

    template<typename T>
    /**
     * @brief Applies the grade involution to a multivector.
     * @param m The multivector to involve.
     * @return The multivector with its odd grades, 1 and 3, negated.
     */
    constexpr Multivector3<T> involute(const Multivector3<T>& m) {
        return Multivector3<T>(
            m.scalar,
            Vector3<T>(-m.vector.x, -m.vector.y, -m.vector.z),
            m.bivector,
            Trivector3<T>(-m.trivector.e123)
        );
    }

    template<typename T>
    /**
     * @brief Applies the grade involution to a rotor.
     * @param r The rotor to involve.
     * @return The rotor unchanged. A rotor holds only even grades, which is
     *         exactly what makes rotors closed under the sandwich product.
     */
    constexpr Rotor3<T> involute(const Rotor3<T>& r) {
        return r;
    }

    // -----------------------------------------------------------------------
    // Clifford conjugation: reverse composed with the grade involution.
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Applies Clifford conjugation to a scalar.
     * @param s The scalar to conjugate.
     * @return The scalar unchanged.
     */
    constexpr Scalar<T> conjugate(const Scalar<T>& s) {
        return s;
    }

    template<typename T>
    /**
     * @brief Applies Clifford conjugation to a vector.
     * @param v The vector to conjugate.
     * @return The negated vector.
     */
    constexpr Vector3<T> conjugate(const Vector3<T>& v) {
        return Vector3<T>(-v.x, -v.y, -v.z);
    }

    template<typename T>
    /**
     * @brief Applies Clifford conjugation to a bivector.
     * @param b The bivector to conjugate.
     * @return The negated bivector.
     */
    constexpr Bivector3<T> conjugate(const Bivector3<T>& b) {
        return Bivector3<T>(-b.xy, -b.xz, -b.yz);
    }

    template<typename T>
    /**
     * @brief Applies Clifford conjugation to a trivector.
     * @param t The trivector to conjugate.
     * @return The trivector unchanged; grades 0 and 3 survive conjugation.
     */
    constexpr Trivector3<T> conjugate(const Trivector3<T>& t) {
        return t;
    }

    template<typename T>
    /**
     * @brief Applies Clifford conjugation to a multivector.
     * @param m The multivector to conjugate.
     * @return The multivector with grades 1 and 2 negated.
     *
     * This is the involution that makes the general inverse work: m multiplied
     * by its conjugate collapses to a scalar plus a pseudoscalar, which then
     * inverts like a complex number. See inverse.hpp.
     */
    constexpr Multivector3<T> conjugate(const Multivector3<T>& m) {
        return Multivector3<T>(
            m.scalar,
            Vector3<T>(-m.vector.x, -m.vector.y, -m.vector.z),
            Bivector3<T>(-m.bivector.xy, -m.bivector.xz, -m.bivector.yz),
            m.trivector
        );
    }

    template<typename T>
    /**
     * @brief Applies Clifford conjugation to a rotor.
     * @param r The rotor to conjugate.
     * @return The rotor with its bivector part negated. For a rotor this
     *         coincides with reverse(), since a rotor has no odd grades.
     */
    constexpr Rotor3<T> conjugate(const Rotor3<T>& r) {
        return Rotor3<T>(r.scalar, Bivector3<T>(-r.bivector.xy, -r.bivector.xz, -r.bivector.yz));
    }
} // namespace CliffordCore
