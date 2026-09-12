#pragma once

/**
 * @file cliffordcore/cl3/operations/reverse.hpp
 * @brief The reverse involution, which flips grades 2 and 3.
 *
 * Reverse scales grade k by (-1)^(k(k-1)/2): grades 0 and 1 keep their
 * sign, grades 2 and 3 are negated. It is the involution rotors are undone
 * with, since inverse(r) is reverse(r) for a unit rotor.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

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
     * @brief Computes the geometric reverse of a scalar.
     * @param s The scalar for which to compute the reverse.
     * @return The scalar unchanged. Grade 0 keeps its sign under reverse.
     */
    constexpr Scalar<T> reverse(const Scalar<T>& s) {
        return s;
    }

    template<typename T>
    /**
     * @brief Computes the geometric reverse of a 3D vector.
     * @param v The vector for which to compute the reverse.
     * @return The reverse of the vector v. (reverse(v) = v)
     */
    constexpr Vector<T> reverse(const Vector<T>& v) {
        return Vector<T>(v.x, v.y, v.z);
    }

    template<typename T>
    /**
     * @brief Computes the geometric reverse of a bivector in 3D space.
     * @param b The bivector for which to compute the reverse.
     * @return The reverse of the bivector b. (reverse(b) = -b)
     */
    constexpr Bivector<T> reverse(const Bivector<T>& b) {
        return Bivector<T>(-b.xy, -b.xz, -b.yz);
    }

    template<typename T>
    /**
     * @brief Computes the geometric reverse of a trivector in 3D space.
     * @param t The trivector for which to compute the reverse.
     * @return The reverse of the trivector t. (reverse(t) = -t)
     */
    constexpr Trivector<T> reverse(const Trivector<T>& t) {
        return Trivector<T>(-t.e123);
    }

    template<typename T>
    /**
     * @brief Computes the geometric reverse of a multivector in 3D space.
     * @param m The multivector for which to compute the reverse.
     * @return The reverse of the multivector m. (reverse(m) = scalar(m) + vector(m) - bivector(m) - trivector(m))
     */
    constexpr Multivector<T> reverse(const Multivector<T>& m) {
        return Multivector<T>(
            m.scalar,
            m.vector,
            Bivector<T>(-m.bivector.xy, -m.bivector.xz, -m.bivector.yz),
            Trivector<T>(-m.trivector.e123)
        );
    }

    template<typename T>
    /**
     * @brief Computes the geometric reverse of a rotor in 3D space.
     * @param r The rotor for which to compute the reverse.
     * @return The reverse of the rotor r. (reverse(r) = scalar(r) - bivector(r))
     */
    constexpr Rotor<T> reverse(const Rotor<T>& r) {
        return Rotor<T>(r.scalar, Bivector<T>(-r.bivector.xy, -r.bivector.xz, -r.bivector.yz));
    }

    


} // namespace CliffordCore::Cl3