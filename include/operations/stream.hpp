#pragma once

/**
 * @file stream.hpp
 * @brief Stream insertion for every type, delegating to to_string().
 */

#include <ostream>

#include "../scalar.hpp"
#include "../vector3.hpp"
#include "../bivector3.hpp"
#include "../trivector3.hpp"
#include "../multivector3.hpp"
#include "../rotor3.hpp"

// Kept out of the type headers so they do not all drag in <ostream>. Include
// this one, or the umbrella header, when you want to print.

namespace CliffordCore
{
    template<typename T>
    /**
     * @brief Writes a scalar to a stream.
     * @param os The stream to write to.
     * @param s The scalar to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Scalar<T>& s) {
        return os << s.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a vector to a stream.
     * @param os The stream to write to.
     * @param v The vector to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Vector3<T>& v) {
        return os << v.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a bivector to a stream.
     * @param os The stream to write to.
     * @param b The bivector to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Bivector3<T>& b) {
        return os << b.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a trivector to a stream.
     * @param os The stream to write to.
     * @param t The trivector to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Trivector3<T>& t) {
        return os << t.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a multivector to a stream.
     * @param os The stream to write to.
     * @param m The multivector to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Multivector3<T>& m) {
        return os << m.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a rotor to a stream.
     * @param os The stream to write to.
     * @param r The rotor to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Rotor3<T>& r) {
        return os << r.to_string();
    }
} // namespace CliffordCore
