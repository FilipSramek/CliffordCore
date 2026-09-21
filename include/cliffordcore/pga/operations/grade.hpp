#pragma once

/**
 * @file cliffordcore/pga/operations/grade.hpp
 * @brief Grade projection: pulling one grade out of a mixed object.
 *
 * grade0 through grade4 pull a single grade out of a mixed object. They are
 * how the closed-form operations extract their result after routing through
 * the general product, and how a rotor or translator hands its bivector half
 * back as a Bivector.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../quadvector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"
#include "../translator.hpp"
#include "../motor.hpp"

// Grade projection: pull a single grade out of a mixed object.
//
// These are named per grade rather than written as grade<N>(m) because each
// grade has a different return type, which a single template cannot express
// without extra machinery.

namespace CliffordCore::PGA
{
    template<typename T>
    /**
     * @brief Extracts the grade 0 (scalar) part of a multivector.
     * @param m The multivector to project.
     * @return The scalar part.
     */
    constexpr Scalar<T> grade0(const Multivector<T>& m) {
        return m.scalar;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 1 (vector, plane) part of a multivector.
     * @param m The multivector to project.
     * @return The vector part.
     */
    constexpr Vector<T> grade1(const Multivector<T>& m) {
        return m.vector;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 2 (bivector, line) part of a multivector.
     * @param m The multivector to project.
     * @return The bivector part.
     */
    constexpr Bivector<T> grade2(const Multivector<T>& m) {
        return m.bivector;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 3 (trivector, point) part of a multivector.
     * @param m The multivector to project.
     * @return The trivector part.
     */
    constexpr Trivector<T> grade3(const Multivector<T>& m) {
        return m.trivector;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 4 (quadvector, pseudoscalar) part of a multivector.
     * @param m The multivector to project.
     * @return The quadvector part.
     */
    constexpr Quadvector<T> grade4(const Multivector<T>& m) {
        return m.quadvector;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 0 (scalar) part of a rotor.
     * @param r The rotor to project.
     * @return The scalar part.
     */
    constexpr Scalar<T> grade0(const Rotor<T>& r) {
        return r.scalar;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 2 (bivector) part of a rotor.
     * @param r The rotor to project.
     * @return The bivector part, with zero ideal components.
     */
    constexpr Bivector<T> grade2(const Rotor<T>& r) {
        return Bivector<T>(0, 0, 0, r.e12, r.e13, r.e23);
    }

    template<typename T>
    /**
     * @brief Extracts the grade 0 (scalar) part of a translator.
     * @param t The translator to project.
     * @return The scalar part.
     */
    constexpr Scalar<T> grade0(const Translator<T>& t) {
        return t.scalar;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 2 (bivector) part of a translator.
     * @param t The translator to project.
     * @return The bivector part, with zero Euclidean components.
     */
    constexpr Bivector<T> grade2(const Translator<T>& t) {
        return Bivector<T>(t.e01, t.e02, t.e03, 0, 0, 0);
    }

    template<typename T>
    /**
     * @brief Extracts the grade 0 (scalar) part of a motor.
     * @param m The motor to project.
     * @return The scalar part.
     */
    constexpr Scalar<T> grade0(const Motor<T>& m) {
        return m.scalar;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 2 (bivector) part of a motor.
     * @param m The motor to project.
     * @return The bivector part.
     */
    constexpr Bivector<T> grade2(const Motor<T>& m) {
        return m.bivector;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 4 (quadvector) part of a motor.
     * @param m The motor to project.
     * @return The quadvector part.
     */
    constexpr Quadvector<T> grade4(const Motor<T>& m) {
        return m.quadvector;
    }
} // namespace CliffordCore::PGA
