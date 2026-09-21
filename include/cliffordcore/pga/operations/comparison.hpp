#pragma once

/**
 * @file cliffordcore/pga/operations/comparison.hpp
 * @brief Exact equality, and a separate tolerant comparison.
 *
 * operator== is exact and usable in a constant expression. approx_equal is
 * kept separate on purpose: a tolerant equality is not transitive, which
 * breaks sorting and associative containers. Both compare component by
 * component, so a motor and its negation -- which represent the same
 * motion -- compare unequal; that is the same choice Cl3 makes for rotors.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <cmath>
#include <limits>

#include "../../detail/compare.hpp"
#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../quadvector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"
#include "../translator.hpp"
#include "../motor.hpp"

// operator== is EXACT, component by component. approx_equal() is separate and
// takes a tolerance. See cl3/operations/comparison.hpp for why the split is
// deliberate: a tolerant == is not an equivalence relation.

namespace CliffordCore::PGA
{
    template<typename T>
    /**
     * @brief Compares two scalars exactly.
     * @param a The left scalar.
     * @param b The right scalar.
     * @return True if every component is bit-identical.
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
        return a.e0 == b.e0
            && a.e1 == b.e1
            && a.e2 == b.e2
            && a.e3 == b.e3;
    }

    template<typename T>
    /**
     * @brief Compares two bivectors exactly.
     * @param a The left bivector.
     * @param b The right bivector.
     * @return True if every component is bit-identical.
     */
    constexpr bool operator==(const Bivector<T>& a, const Bivector<T>& b) {
        return a.e01 == b.e01
            && a.e02 == b.e02
            && a.e03 == b.e03
            && a.e12 == b.e12
            && a.e13 == b.e13
            && a.e23 == b.e23;
    }

    template<typename T>
    /**
     * @brief Compares two trivectors exactly.
     * @param a The left trivector.
     * @param b The right trivector.
     * @return True if every component is bit-identical.
     */
    constexpr bool operator==(const Trivector<T>& a, const Trivector<T>& b) {
        return a.e012 == b.e012
            && a.e013 == b.e013
            && a.e023 == b.e023
            && a.e123 == b.e123;
    }

    template<typename T>
    /**
     * @brief Compares two quadvectors exactly.
     * @param a The left quadvector.
     * @param b The right quadvector.
     * @return True if every component is bit-identical.
     */
    constexpr bool operator==(const Quadvector<T>& a, const Quadvector<T>& b) {
        return a.e0123 == b.e0123;
    }

    template<typename T>
    /**
     * @brief Compares two multivectors exactly.
     * @param a The left multivector.
     * @param b The right multivector.
     * @return True if every component is bit-identical.
     */
    constexpr bool operator==(const Multivector<T>& a, const Multivector<T>& b) {
        return a.scalar == b.scalar
            && a.vector == b.vector
            && a.bivector == b.bivector
            && a.trivector == b.trivector
            && a.quadvector == b.quadvector;
    }

    template<typename T>
    /**
     * @brief Compares two rotors exactly.
     * @param a The left rotor.
     * @param b The right rotor.
     * @return True if every component is bit-identical.
     */
    constexpr bool operator==(const Rotor<T>& a, const Rotor<T>& b) {
        return a.scalar.value == b.scalar.value
            && a.e12 == b.e12
            && a.e13 == b.e13
            && a.e23 == b.e23;
    }

    template<typename T>
    /**
     * @brief Compares two translators exactly.
     * @param a The left translator.
     * @param b The right translator.
     * @return True if every component is bit-identical.
     */
    constexpr bool operator==(const Translator<T>& a, const Translator<T>& b) {
        return a.scalar.value == b.scalar.value
            && a.e01 == b.e01
            && a.e02 == b.e02
            && a.e03 == b.e03;
    }

    template<typename T>
    /**
     * @brief Compares two motors exactly.
     * @param a The left motor.
     * @param b The right motor.
     * @return True if every component is bit-identical.
     */
    constexpr bool operator==(const Motor<T>& a, const Motor<T>& b) {
        return a.scalar == b.scalar
            && a.bivector == b.bivector
            && a.quadvector == b.quadvector;
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
    constexpr bool operator!=(const Trivector<T>& a, const Trivector<T>& b) { return !(a == b); }

    template<typename T>
    /** @brief Negation of operator==. @param a Left. @param b Right. @return True if not exactly equal. */
    constexpr bool operator!=(const Quadvector<T>& a, const Quadvector<T>& b) { return !(a == b); }

    template<typename T>
    /** @brief Negation of operator==. @param a Left. @param b Right. @return True if not exactly equal. */
    constexpr bool operator!=(const Multivector<T>& a, const Multivector<T>& b) { return !(a == b); }

    template<typename T>
    /** @brief Negation of operator==. @param a Left. @param b Right. @return True if not exactly equal. */
    constexpr bool operator!=(const Rotor<T>& a, const Rotor<T>& b) { return !(a == b); }

    template<typename T>
    /** @brief Negation of operator==. @param a Left. @param b Right. @return True if not exactly equal. */
    constexpr bool operator!=(const Translator<T>& a, const Translator<T>& b) { return !(a == b); }

    template<typename T>
    /** @brief Negation of operator==. @param a Left. @param b Right. @return True if not exactly equal. */
    constexpr bool operator!=(const Motor<T>& a, const Motor<T>& b) { return !(a == b); }

    template<typename T>
    /**
     * @brief Compares two scalars within a tolerance.
     * @param a The left scalar.
     * @param b The right scalar.
     * @param tolerance The largest per-component difference still considered equal.
     * @return True if every component agrees to within the tolerance.
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
        return CliffordCore::detail::close(a.e0, b.e0, tolerance)
            && CliffordCore::detail::close(a.e1, b.e1, tolerance)
            && CliffordCore::detail::close(a.e2, b.e2, tolerance)
            && CliffordCore::detail::close(a.e3, b.e3, tolerance);
    }

    template<typename T>
    /**
     * @brief Compares two bivectors within a tolerance.
     * @param a The left bivector.
     * @param b The right bivector.
     * @param tolerance The largest per-component difference still considered equal.
     * @return True if every component agrees to within the tolerance.
     */
    inline bool approx_equal(const Bivector<T>& a, const Bivector<T>& b,
                             T tolerance = CliffordCore::detail::default_tolerance<T>()) {
        return CliffordCore::detail::close(a.e01, b.e01, tolerance)
            && CliffordCore::detail::close(a.e02, b.e02, tolerance)
            && CliffordCore::detail::close(a.e03, b.e03, tolerance)
            && CliffordCore::detail::close(a.e12, b.e12, tolerance)
            && CliffordCore::detail::close(a.e13, b.e13, tolerance)
            && CliffordCore::detail::close(a.e23, b.e23, tolerance);
    }

    template<typename T>
    /**
     * @brief Compares two trivectors within a tolerance.
     * @param a The left trivector.
     * @param b The right trivector.
     * @param tolerance The largest per-component difference still considered equal.
     * @return True if every component agrees to within the tolerance.
     */
    inline bool approx_equal(const Trivector<T>& a, const Trivector<T>& b,
                             T tolerance = CliffordCore::detail::default_tolerance<T>()) {
        return CliffordCore::detail::close(a.e012, b.e012, tolerance)
            && CliffordCore::detail::close(a.e013, b.e013, tolerance)
            && CliffordCore::detail::close(a.e023, b.e023, tolerance)
            && CliffordCore::detail::close(a.e123, b.e123, tolerance);
    }

    template<typename T>
    /**
     * @brief Compares two quadvectors within a tolerance.
     * @param a The left quadvector.
     * @param b The right quadvector.
     * @param tolerance The largest per-component difference still considered equal.
     * @return True if every component agrees to within the tolerance.
     */
    inline bool approx_equal(const Quadvector<T>& a, const Quadvector<T>& b,
                             T tolerance = CliffordCore::detail::default_tolerance<T>()) {
        return CliffordCore::detail::close(a.e0123, b.e0123, tolerance);
    }

    template<typename T>
    /**
     * @brief Compares two multivectors within a tolerance.
     * @param a The left multivector.
     * @param b The right multivector.
     * @param tolerance The largest per-component difference still considered equal.
     * @return True if every component agrees to within the tolerance.
     */
    inline bool approx_equal(const Multivector<T>& a, const Multivector<T>& b,
                             T tolerance = CliffordCore::detail::default_tolerance<T>()) {
        return approx_equal(a.scalar, b.scalar, tolerance)
            && approx_equal(a.vector, b.vector, tolerance)
            && approx_equal(a.bivector, b.bivector, tolerance)
            && approx_equal(a.trivector, b.trivector, tolerance)
            && approx_equal(a.quadvector, b.quadvector, tolerance);
    }

    template<typename T>
    /**
     * @brief Compares two rotors within a tolerance.
     * @param a The left rotor.
     * @param b The right rotor.
     * @param tolerance The largest per-component difference still considered equal.
     * @return True if every component agrees to within the tolerance.
     */
    inline bool approx_equal(const Rotor<T>& a, const Rotor<T>& b,
                             T tolerance = CliffordCore::detail::default_tolerance<T>()) {
        return CliffordCore::detail::close(a.scalar.value, b.scalar.value, tolerance)
            && CliffordCore::detail::close(a.e12, b.e12, tolerance)
            && CliffordCore::detail::close(a.e13, b.e13, tolerance)
            && CliffordCore::detail::close(a.e23, b.e23, tolerance);
    }

    template<typename T>
    /**
     * @brief Compares two translators within a tolerance.
     * @param a The left translator.
     * @param b The right translator.
     * @param tolerance The largest per-component difference still considered equal.
     * @return True if every component agrees to within the tolerance.
     */
    inline bool approx_equal(const Translator<T>& a, const Translator<T>& b,
                             T tolerance = CliffordCore::detail::default_tolerance<T>()) {
        return CliffordCore::detail::close(a.scalar.value, b.scalar.value, tolerance)
            && CliffordCore::detail::close(a.e01, b.e01, tolerance)
            && CliffordCore::detail::close(a.e02, b.e02, tolerance)
            && CliffordCore::detail::close(a.e03, b.e03, tolerance);
    }

    template<typename T>
    /**
     * @brief Compares two motors within a tolerance.
     * @param a The left motor.
     * @param b The right motor.
     * @param tolerance The largest per-component difference still considered equal.
     * @return True if every component agrees to within the tolerance.
     */
    inline bool approx_equal(const Motor<T>& a, const Motor<T>& b,
                             T tolerance = CliffordCore::detail::default_tolerance<T>()) {
        return approx_equal(a.scalar, b.scalar, tolerance)
            && approx_equal(a.bivector, b.bivector, tolerance)
            && approx_equal(a.quadvector, b.quadvector, tolerance);
    }
} // namespace CliffordCore::PGA
