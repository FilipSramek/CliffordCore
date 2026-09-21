#pragma once

/**
 * @file cliffordcore/pga/operations/sandwich.hpp
 * @brief Applying a motor to geometry via the sandwich product.
 *
 * sandwich(x, m) applies M x ~M, with the object first and the versor
 * second -- the same argument order as Cl3's. One formula moves planes,
 * lines, points and whole multivectors, which is the point of PGA: a
 * rigid motion acts on every kind of object the same way.
 *
 * rotate, translate and transform are the same function under friendlier
 * names for the three versor types. The Rotor and Translator overloads
 * widen to a Motor first, so there is one implementation to get right.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <type_traits>
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"
#include "../translator.hpp"
#include "../motor.hpp"
#include "grade.hpp"
#include "mixed_products.hpp"
#include "reverse.hpp"

namespace CliffordCore::PGA
{
    template<typename T>
    /**
     * @brief Applies a motor to a multivector: M x ~M.
     * @param x The multivector to transform.
     * @param m The motor to transform it with.
     * @return The transformed multivector.
     */
    constexpr Multivector<T> sandwich(const Multivector<T>& x, const Motor<T>& m) {
        return detail::promote(m) * x * detail::promote(reverse(m));
    }

    template<typename T>
    /**
     * @brief Applies a motor to a plane: M p ~M.
     * @param x The vector (plane) to transform.
     * @param m The motor to transform it with.
     * @return The transformed plane. A motor is even, so the grade is preserved.
     */
    constexpr Vector<T> sandwich(const Vector<T>& x, const Motor<T>& m) {
        return grade1(sandwich(detail::promote(x), m));
    }

    template<typename T>
    /**
     * @brief Applies a motor to a line: M L ~M.
     * @param x The bivector (line) to transform.
     * @param m The motor to transform it with.
     * @return The transformed line.
     */
    constexpr Bivector<T> sandwich(const Bivector<T>& x, const Motor<T>& m) {
        return grade2(sandwich(detail::promote(x), m));
    }

    template<typename T>
    /**
     * @brief Applies a motor to a point: M P ~M.
     * @param x The trivector (point) to transform.
     * @param m The motor to transform it with.
     * @return The transformed point, with its weight preserved.
     */
    constexpr Trivector<T> sandwich(const Trivector<T>& x, const Motor<T>& m) {
        return grade3(sandwich(detail::promote(x), m));
    }

    // Rotor and Translator overloads. Each widens to a Motor, which is lossless.

    template<typename T>
    /**
     * @brief Applies a rotor to a plane.
     * @param x The vector (plane) to transform.
     * @param r The rotor to transform it with.
     * @return The transformed plane.
     */
    constexpr Vector<T> sandwich(const Vector<T>& x, const Rotor<T>& r) {
        return sandwich(x, Motor<T>(r));
    }

    template<typename T>
    /**
     * @brief Applies a rotor to a line.
     * @param x The bivector (line) to transform.
     * @param r The rotor to transform it with.
     * @return The transformed line.
     */
    constexpr Bivector<T> sandwich(const Bivector<T>& x, const Rotor<T>& r) {
        return sandwich(x, Motor<T>(r));
    }

    template<typename T>
    /**
     * @brief Applies a rotor to a point.
     * @param x The trivector (point) to transform.
     * @param r The rotor to transform it with.
     * @return The transformed point.
     */
    constexpr Trivector<T> sandwich(const Trivector<T>& x, const Rotor<T>& r) {
        return sandwich(x, Motor<T>(r));
    }

    template<typename T>
    /**
     * @brief Applies a rotor to a multivector.
     * @param x The multivector to transform.
     * @param r The rotor to transform it with.
     * @return The transformed multivector.
     */
    constexpr Multivector<T> sandwich(const Multivector<T>& x, const Rotor<T>& r) {
        return sandwich(x, Motor<T>(r));
    }

    template<typename T>
    /**
     * @brief Applies a translator to a plane.
     * @param x The vector (plane) to transform.
     * @param t The translator to transform it with.
     * @return The translated plane: its normal is unchanged, its offset moves.
     */
    constexpr Vector<T> sandwich(const Vector<T>& x, const Translator<T>& t) {
        return sandwich(x, Motor<T>(t));
    }

    template<typename T>
    /**
     * @brief Applies a translator to a line.
     * @param x The bivector (line) to transform.
     * @param t The translator to transform it with.
     * @return The translated line: its direction is unchanged, its moment moves.
     */
    constexpr Bivector<T> sandwich(const Bivector<T>& x, const Translator<T>& t) {
        return sandwich(x, Motor<T>(t));
    }

    template<typename T>
    /**
     * @brief Applies a translator to a point.
     * @param x The trivector (point) to transform.
     * @param t The translator to transform it with.
     * @return The translated point.
     */
    constexpr Trivector<T> sandwich(const Trivector<T>& x, const Translator<T>& t) {
        return sandwich(x, Motor<T>(t));
    }

    template<typename T>
    /**
     * @brief Applies a translator to a multivector.
     * @param x The multivector to transform.
     * @param t The translator to transform it with.
     * @return The translated multivector.
     */
    constexpr Multivector<T> sandwich(const Multivector<T>& x, const Translator<T>& t) {
        return sandwich(x, Motor<T>(t));
    }

    // -----------------------------------------------------------------------
    // Friendlier names for the three versor types
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Rotates a plane about the rotor's axis through the origin.
     * @param x The vector (plane) to rotate.
     * @param r The rotor.
     * @return sandwich(x, r).
     */
    constexpr Vector<T> rotate(const Vector<T>& x, const Rotor<T>& r) {
        return sandwich(x, r);
    }

    template<typename T>
    /**
     * @brief Rotates a line about the rotor's axis through the origin.
     * @param x The bivector (line) to rotate.
     * @param r The rotor.
     * @return sandwich(x, r).
     */
    constexpr Bivector<T> rotate(const Bivector<T>& x, const Rotor<T>& r) {
        return sandwich(x, r);
    }

    template<typename T>
    /**
     * @brief Rotates a point about the rotor's axis through the origin.
     * @param x The trivector (point) to rotate.
     * @param r The rotor.
     * @return sandwich(x, r).
     */
    constexpr Trivector<T> rotate(const Trivector<T>& x, const Rotor<T>& r) {
        return sandwich(x, r);
    }

    template<typename T>
    /**
     * @brief Rotates a multivector about the rotor's axis through the origin.
     * @param x The multivector to rotate.
     * @param r The rotor.
     * @return sandwich(x, r).
     */
    constexpr Multivector<T> rotate(const Multivector<T>& x, const Rotor<T>& r) {
        return sandwich(x, r);
    }

    template<typename T>
    /**
     * @brief Translates a plane.
     * @param x The vector (plane) to translate.
     * @param t The translator.
     * @return sandwich(x, t).
     */
    constexpr Vector<T> translate(const Vector<T>& x, const Translator<T>& t) {
        return sandwich(x, t);
    }

    template<typename T>
    /**
     * @brief Translates a line.
     * @param x The bivector (line) to translate.
     * @param t The translator.
     * @return sandwich(x, t).
     */
    constexpr Bivector<T> translate(const Bivector<T>& x, const Translator<T>& t) {
        return sandwich(x, t);
    }

    template<typename T>
    /**
     * @brief Translates a point.
     * @param x The trivector (point) to translate.
     * @param t The translator.
     * @return sandwich(x, t).
     */
    constexpr Trivector<T> translate(const Trivector<T>& x, const Translator<T>& t) {
        return sandwich(x, t);
    }

    template<typename T>
    /**
     * @brief Translates a multivector.
     * @param x The multivector to translate.
     * @param t The translator.
     * @return sandwich(x, t).
     */
    constexpr Multivector<T> translate(const Multivector<T>& x, const Translator<T>& t) {
        return sandwich(x, t);
    }

    template<typename T>
    /**
     * @brief Applies a rigid motion to a plane.
     * @param x The vector (plane) to transform.
     * @param m The motor.
     * @return sandwich(x, m).
     */
    constexpr Vector<T> transform(const Vector<T>& x, const Motor<T>& m) {
        return sandwich(x, m);
    }

    template<typename T>
    /**
     * @brief Applies a rigid motion to a line.
     * @param x The bivector (line) to transform.
     * @param m The motor.
     * @return sandwich(x, m).
     */
    constexpr Bivector<T> transform(const Bivector<T>& x, const Motor<T>& m) {
        return sandwich(x, m);
    }

    template<typename T>
    /**
     * @brief Applies a rigid motion to a point.
     * @param x The trivector (point) to transform.
     * @param m The motor.
     * @return sandwich(x, m).
     */
    constexpr Trivector<T> transform(const Trivector<T>& x, const Motor<T>& m) {
        return sandwich(x, m);
    }

    template<typename T>
    /**
     * @brief Applies a rigid motion to a multivector.
     * @param x The multivector to transform.
     * @param m The motor.
     * @return sandwich(x, m).
     */
    constexpr Multivector<T> transform(const Multivector<T>& x, const Motor<T>& m) {
        return sandwich(x, m);
    }
} // namespace CliffordCore::PGA
