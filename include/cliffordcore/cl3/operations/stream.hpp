#pragma once

/**
 * @file cliffordcore/cl3/operations/stream.hpp
 * @brief Stream insertion for every type, delegating to to_string().
 *
 * One-line operator<< for each type, delegating to its to_string(). Kept in
 * its own header so the type headers do not have to pull in <ostream>.
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <ostream>

#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"

// Kept out of the type headers so they do not all drag in <ostream>. Include
// this one, or the umbrella header, when you want to print.

namespace CliffordCore::Cl3
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
    inline std::ostream& operator<<(std::ostream& os, const Vector<T>& v) {
        return os << v.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a bivector to a stream.
     * @param os The stream to write to.
     * @param b The bivector to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Bivector<T>& b) {
        return os << b.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a trivector to a stream.
     * @param os The stream to write to.
     * @param t The trivector to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Trivector<T>& t) {
        return os << t.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a multivector to a stream.
     * @param os The stream to write to.
     * @param m The multivector to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Multivector<T>& m) {
        return os << m.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a rotor to a stream.
     * @param os The stream to write to.
     * @param r The rotor to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Rotor<T>& r) {
        return os << r.to_string();
    }
} // namespace CliffordCore::Cl3
