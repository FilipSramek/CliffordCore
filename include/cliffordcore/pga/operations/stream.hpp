#pragma once

/**
 * @file cliffordcore/pga/operations/stream.hpp
 * @brief Stream insertion for every type, delegating to to_string().
 *
 * One-line operator<< for each type, delegating to its to_string(). Kept in
 * its own header so the type headers do not have to pull in <ostream>.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <ostream>

#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../quadvector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"
#include "../translator.hpp"
#include "../motor.hpp"

// Kept out of the type headers so they do not all drag in <ostream>. Include
// this one, or the umbrella header, when you want to print.

namespace CliffordCore::PGA
{
    template<typename T>
    /**
     * @brief Writes a scalar to a stream.
     * @param os The stream to write to.
     * @param x The scalar to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Scalar<T>& x) {
        return os << x.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a vector to a stream.
     * @param os The stream to write to.
     * @param x The vector to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Vector<T>& x) {
        return os << x.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a bivector to a stream.
     * @param os The stream to write to.
     * @param x The bivector to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Bivector<T>& x) {
        return os << x.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a trivector to a stream.
     * @param os The stream to write to.
     * @param x The trivector to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Trivector<T>& x) {
        return os << x.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a quadvector to a stream.
     * @param os The stream to write to.
     * @param x The quadvector to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Quadvector<T>& x) {
        return os << x.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a multivector to a stream.
     * @param os The stream to write to.
     * @param x The multivector to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Multivector<T>& x) {
        return os << x.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a rotor to a stream.
     * @param os The stream to write to.
     * @param x The rotor to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Rotor<T>& x) {
        return os << x.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a translator to a stream.
     * @param os The stream to write to.
     * @param x The translator to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Translator<T>& x) {
        return os << x.to_string();
    }

    template<typename T>
    /**
     * @brief Writes a motor to a stream.
     * @param os The stream to write to.
     * @param x The motor to write.
     * @return The stream.
     */
    inline std::ostream& operator<<(std::ostream& os, const Motor<T>& x) {
        return os << x.to_string();
    }
} // namespace CliffordCore::PGA
