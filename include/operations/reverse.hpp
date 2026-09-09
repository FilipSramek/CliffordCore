#pragma once

/**
 * @file reverse.hpp
 * @brief The reverse involution, which flips grades 2 and 3.
 */

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
    constexpr Vector3<T> reverse(const Vector3<T>& v) {
        return Vector3<T>(v.x, v.y, v.z);
    }

    template<typename T>
    /**
     * @brief Computes the geometric reverse of a bivector in 3D space.
     * @param b The bivector for which to compute the reverse.
     * @return The reverse of the bivector b. (reverse(b) = -b)
     */
    constexpr Bivector3<T> reverse(const Bivector3<T>& b) {
        return Bivector3<T>(-b.xy, -b.xz, -b.yz);
    }

    template<typename T>
    /**
     * @brief Computes the geometric reverse of a trivector in 3D space.
     * @param t The trivector for which to compute the reverse.
     * @return The reverse of the trivector t. (reverse(t) = -t)
     */
    constexpr Trivector3<T> reverse(const Trivector3<T>& t) {
        return Trivector3<T>(-t.e123);
    }

    template<typename T>
    /**
     * @brief Computes the geometric reverse of a multivector in 3D space.
     * @param m The multivector for which to compute the reverse.
     * @return The reverse of the multivector m. (reverse(m) = scalar(m) + vector(m) - bivector(m) - trivector(m))
     */
    constexpr Multivector3<T> reverse(const Multivector3<T>& m) {
        return Multivector3<T>(
            m.scalar,
            m.vector,
            Bivector3<T>(-m.bivector.xy, -m.bivector.xz, -m.bivector.yz),
            Trivector3<T>(-m.trivector.e123)
        );
    }

    template<typename T>
    /**
     * @brief Computes the geometric reverse of a rotor in 3D space.
     * @param r The rotor for which to compute the reverse.
     * @return The reverse of the rotor r. (reverse(r) = scalar(r) - bivector(r))
     */
    constexpr Rotor3<T> reverse(const Rotor3<T>& r) {
        return Rotor3<T>(r.scalar, Bivector3<T>(-r.bivector.xy, -r.bivector.xz, -r.bivector.yz));
    }

    


} // namespace CliffordCore