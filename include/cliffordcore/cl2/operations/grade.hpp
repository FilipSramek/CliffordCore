#pragma once

/**
 * @file cliffordcore/cl2/operations/grade.hpp
 * @brief Grade projection: pulling one grade out of a mixed object.
 *
 * grade0 through grade2 pull a single grade out of a mixed object. They are
 * how the closed-form operations extract their result after routing through
 * the general product.
 *
 * @author Filip Sramek
 * @version 0.3.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../multivector.hpp"
#include "../rotor.hpp"

// Grade projection: pull a single grade out of a mixed object.
//
// These are named per grade rather than written as grade<N>(m) because each
// grade has a different return type, which a single template cannot express
// without extra machinery.
//
// There is one overload Cl(3,0) does not have: grade0 of a Scalar, which is
// the identity. In Cl(2,0) two bivectors multiply to a pure Scalar, and
// scalar_product(a, b) is spelled grade0(a * b) -- without this overload the
// scalar product of two bivectors would fail to compile.

namespace CliffordCore::Cl2
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
     * @brief Extracts the grade 1 (vector) part of a multivector.
     * @param m The multivector to project.
     * @return The vector part.
     */
    constexpr Vector<T> grade1(const Multivector<T>& m) {
        return m.vector;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 2 (bivector) part of a multivector.
     * @param m The multivector to project.
     * @return The bivector part, which is also the pseudoscalar part.
     */
    constexpr Bivector<T> grade2(const Multivector<T>& m) {
        return m.bivector;
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
     * @return The bivector part.
     */
    constexpr Bivector<T> grade2(const Rotor<T>& r) {
        return r.bivector;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 0 part of a scalar, which is the scalar itself.
     * @param s The scalar to project.
     * @return The same scalar.
     */
    constexpr Scalar<T> grade0(const Scalar<T>& s) {
        return s;
    }
} // namespace CliffordCore::Cl2
