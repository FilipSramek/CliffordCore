#pragma once

#include "../scalar.hpp"
#include "../vector3.hpp"
#include "../bivector3.hpp"
#include "../trivector3.hpp"
#include "../multivector3.hpp"
#include "../rotor3.hpp"

// Grade projection: pull a single grade out of a mixed object.
//
// These are named per grade rather than written as grade<N>(m) because each
// grade has a different return type, which a single template cannot express
// without extra machinery.

namespace CliffordCore
{
    template<typename T>
    /**
     * @brief Extracts the grade 0 (scalar) part of a multivector.
     * @param m The multivector to project.
     * @return The scalar part.
     */
    constexpr Scalar<T> grade0(const Multivector3<T>& m) {
        return m.scalar;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 1 (vector) part of a multivector.
     * @param m The multivector to project.
     * @return The vector part.
     */
    constexpr Vector3<T> grade1(const Multivector3<T>& m) {
        return m.vector;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 2 (bivector) part of a multivector.
     * @param m The multivector to project.
     * @return The bivector part.
     */
    constexpr Bivector3<T> grade2(const Multivector3<T>& m) {
        return m.bivector;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 3 (trivector) part of a multivector.
     * @param m The multivector to project.
     * @return The trivector part.
     */
    constexpr Trivector3<T> grade3(const Multivector3<T>& m) {
        return m.trivector;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 0 (scalar) part of a rotor.
     * @param r The rotor to project.
     * @return The scalar part.
     */
    constexpr Scalar<T> grade0(const Rotor3<T>& r) {
        return r.scalar;
    }

    template<typename T>
    /**
     * @brief Extracts the grade 2 (bivector) part of a rotor.
     * @param r The rotor to project.
     * @return The bivector part.
     */
    constexpr Bivector3<T> grade2(const Rotor3<T>& r) {
        return r.bivector;
    }
} // namespace CliffordCore
