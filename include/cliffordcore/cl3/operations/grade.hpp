#pragma once

/**
 * @file cliffordcore/cl3/operations/grade.hpp
 * @brief Grade projection: pulling one grade out of a mixed object.
 *
 * grade0 through grade3 pull a single grade out of a mixed object. They are
 * how the closed-form operations extract their result after routing through
 * the general product.
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

// Grade projection: pull a single grade out of a mixed object.
//
// These are named per grade rather than written as grade<N>(m) because each
// grade has a different return type, which a single template cannot express
// without extra machinery.

namespace CliffordCore::Cl3
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
     * @return The bivector part.
     */
    constexpr Bivector<T> grade2(const Multivector<T>& m) {
        return m.bivector;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 3 (trivector) part of a multivector.
     * @param m The multivector to project.
     * @return The trivector part.
     */
    constexpr Trivector<T> grade3(const Multivector<T>& m) {
        return m.trivector;
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
} // namespace CliffordCore::Cl3
