#pragma once

/**
 * @file cliffordcore/cl2/operations/comparison.hpp
 * @brief Exact equality, and a separate tolerant comparison.
 *
 * operator== is exact and usable in a constant expression. approx_equal is
 * kept separate on purpose: a tolerant equality is not transitive, which
 * breaks sorting and associative containers.
 *
 * @author Filip Sramek
 * @version 0.3.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <cmath>
#include <limits>

#include "../../detail/compare.hpp"
#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"

// operator== is EXACT, component by component. approx_equal() is separate and
// takes a tolerance.
//
// The split is deliberate. A tolerant operator== is not an equivalence
// relation: with a fixed tolerance you can have a == b and b == c but a != c,
// which quietly breaks sorting, associative containers, and every algorithm
// that assumes transitivity. Exact == keeps those guarantees; when you want
// "close enough", you have to say how close, which is the question you should
// be answering anyway.
//
// Exact == is still useful here: results that are equal by construction really
// are bit-identical. Anything that has been through a rotation or an inverse
// wants approx_equal.

namespace CliffordCore::Cl2
{
    template<typename T>
    /**
     * @brief Compares two scalars exactly.
     * @param a The left scalar.
     * @param b The right scalar.
     * @return True if the values are bit-identical.
     */
    constexpr bool operator==(const Scalar<T>& a, const Scalar<T>& b) {
        return a.value == b.value;
    }

    template<typename T>
    /**
     * @brief Compares two vectors exactly.
     * @param a The left vector.
     * @param b The right vector.
     * @return True if every component is bit-identical.
     */
    constexpr bool operator==(const Vector<T>& a, const Vector<T>& b) {
        return a.x == b.x && a.y == b.y;
    }

    template<typename T>
    /**
     * @brief Compares two bivectors exactly.
     * @param a The left bivector.
     * @param b The right bivector.
     * @return True if the components are bit-identical.
     */
    constexpr bool operator==(const Bivector<T>& a, const Bivector<T>& b) {
        return a.xy == b.xy;
    }

    template<typename T>
    /**
     * @brief Compares two multivectors exactly.
     * @param a The left multivector.
     * @param b The right multivector.
     * @return True if every grade is bit-identical.
     */
    constexpr bool operator==(const Multivector<T>& a, const Multivector<T>& b) {
        return a.scalar == b.scalar && a.vector == b.vector && a.bivector == b.bivector;
    }

    template<typename T>
    /**
     * @brief Compares two rotors exactly.
     * @param a The left rotor.
     * @param b The right rotor.
     * @return True if both parts are bit-identical. Note that a rotor and its
     *         negation describe the same rotation but are NOT equal here.
     */
    constexpr bool operator==(const Rotor<T>& a, const Rotor<T>& b) {
        return a.scalar == b.scalar && a.bivector == b.bivector;
    }

    template<typename T>
    /** @brief Negation of operator==. @param a Left. @param b Right. @return True if not exactly equal. */
    constexpr bool operator!=(const Scalar<T>& a, const Scalar<T>& b) { return !(a == b); }

    template<typename T>
    /** @brief Negation of operator==. @param a Left. @param b Right. @return True if not exactly equal. */
    constexpr bool operator!=(const Vector<T>& a, const Vector<T>& b) { return !(a == b); }

    template<typename T>
    /** @brief Negation of operator==. @param a Left. @param b Right. @return True if not exactly equal. */
    constexpr bool operator!=(const Bivector<T>& a, const Bivector<T>& b) { return !(a == b); }

    template<typename T>
    /** @brief Negation of operator==. @param a Left. @param b Right. @return True if not exactly equal. */
    constexpr bool operator!=(const Multivector<T>& a, const Multivector<T>& b) { return !(a == b); }

    template<typename T>
    /** @brief Negation of operator==. @param a Left. @param b Right. @return True if not exactly equal. */
    constexpr bool operator!=(const Rotor<T>& a, const Rotor<T>& b) { return !(a == b); }

    // -----------------------------------------------------------------------
    // Tolerant comparison
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Compares two scalars within a tolerance.
     * @param a The left scalar.
     * @param b The right scalar.
     * @param tolerance The largest difference still considered equal.
     * @return True if they agree to within the tolerance.
     */
    inline bool approx_equal(const Scalar<T>& a, const Scalar<T>& b,
                             T tolerance = CliffordCore::detail::default_tolerance<T>()) {
        return CliffordCore::detail::close(a.value, b.value, tolerance);
    }

    template<typename T>
    /**
     * @brief Compares two vectors within a tolerance.
     * @param a The left vector.
     * @param b The right vector.
     * @param tolerance The largest per-component difference still considered equal.
     * @return True if every component agrees to within the tolerance.
     */
    inline bool approx_equal(const Vector<T>& a, const Vector<T>& b,
                             T tolerance = CliffordCore::detail::default_tolerance<T>()) {
        return CliffordCore::detail::close(a.x, b.x, tolerance)
            && CliffordCore::detail::close(a.y, b.y, tolerance);
    }

    template<typename T>
    /**
     * @brief Compares two bivectors within a tolerance.
     * @param a The left bivector.
     * @param b The right bivector.
     * @param tolerance The largest difference still considered equal.
     * @return True if they agree to within the tolerance.
     */
    inline bool approx_equal(const Bivector<T>& a, const Bivector<T>& b,
                             T tolerance = CliffordCore::detail::default_tolerance<T>()) {
        return CliffordCore::detail::close(a.xy, b.xy, tolerance);
    }

    template<typename T>
    /**
     * @brief Compares two multivectors within a tolerance.
     * @param a The left multivector.
     * @param b The right multivector.
     * @param tolerance The largest per-component difference still considered equal.
     * @return True if every grade agrees to within the tolerance.
     */
    inline bool approx_equal(const Multivector<T>& a, const Multivector<T>& b,
                             T tolerance = CliffordCore::detail::default_tolerance<T>()) {
        return approx_equal(a.scalar, b.scalar, tolerance)
            && approx_equal(a.vector, b.vector, tolerance)
            && approx_equal(a.bivector, b.bivector, tolerance);
    }

    template<typename T>
    /**
     * @brief Compares two rotors within a tolerance.
     * @param a The left rotor.
     * @param b The right rotor.
     * @param tolerance The largest per-component difference still considered equal.
     * @return True if both parts agree to within the tolerance. A rotor and its
     *         negation are the same rotation but are not approx_equal.
     */
    inline bool approx_equal(const Rotor<T>& a, const Rotor<T>& b,
                             T tolerance = CliffordCore::detail::default_tolerance<T>()) {
        return approx_equal(a.scalar, b.scalar, tolerance)
            && approx_equal(a.bivector, b.bivector, tolerance);
    }
} // namespace CliffordCore::Cl2
