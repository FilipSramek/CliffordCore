#pragma once

/**
 * @file cliffordcore/pga/multivector.hpp
 * @brief All grades at once: the Multivector type of Cl(3,0,1).
 *
 * Multivector<T> holds all five grades -- sixteen components -- and is the
 * return type of any product that can leave its operands' grades. Mixed-
 * grade addition and subtraction always widen to it, so chaining stays
 * predictable. Rotor, Translator and Motor widen to it implicitly, because
 * nothing is lost.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <type_traits>
#include <string>
#include "scalar.hpp"
#include "vector.hpp"
#include "bivector.hpp"
#include "trivector.hpp"
#include "quadvector.hpp"
#include "rotor.hpp"
#include "translator.hpp"
#include "motor.hpp"

namespace CliffordCore::PGA
{
    template<typename T>
    /**
     * @brief A class representing a general element of Cl(3,0,1).
     * @tparam T The arithmetic component type.
     */
    class Multivector
    {
        static_assert(std::is_arithmetic<T>::value, "Multivector can only be instantiated with numeric types.");
    public:
        Scalar<T> scalar;          ///< The grade 0 part.
        Vector<T> vector;          ///< The grade 1 part.
        Bivector<T> bivector;      ///< The grade 2 part.
        Trivector<T> trivector;    ///< The grade 3 part.
        Quadvector<T> quadvector;  ///< The grade 4 part.

        /**
         * @brief Default constructor initializes all components to zero.
         */
        constexpr Multivector() : scalar(), vector(), bivector(), trivector(), quadvector() {}

        /**
         * @brief Constructor initializes the multivector with the provided components.
         * @param s The scalar component.
         * @param v The vector component.
         * @param b The bivector component.
         * @param t The trivector component.
         * @param q The quadvector component.
         */
        constexpr Multivector(const Scalar<T>& s, const Vector<T>& v, const Bivector<T>& b,
                              const Trivector<T>& t, const Quadvector<T>& q)
            : scalar(s), vector(v), bivector(b), trivector(t), quadvector(q) {}

        /**
         * @brief Constructor widens a rotor to a multivector.
         * @param r The rotor to widen.
         */
        constexpr Multivector(const Rotor<T>& r)
            : scalar(r.scalar), vector(), bivector(Bivector<T>(0, 0, 0, r.e12, r.e13, r.e23)),
              trivector(), quadvector() {}

        /**
         * @brief Constructor widens a translator to a multivector.
         * @param t The translator to widen.
         */
        constexpr Multivector(const Translator<T>& t)
            : scalar(t.scalar), vector(), bivector(Bivector<T>(t.e01, t.e02, t.e03, 0, 0, 0)),
              trivector(), quadvector() {}

        /**
         * @brief Constructor widens a motor to a multivector.
         * @param m The motor to widen.
         */
        constexpr Multivector(const Motor<T>& m)
            : scalar(m.scalar), vector(), bivector(m.bivector), trivector(), quadvector(m.quadvector) {}

        /**
         * @brief Addition operator overload.
         * @param other The multivector to add.
         * @return The resulting multivector.
         */
        constexpr Multivector operator+(const Multivector& other) const {
            return Multivector(
                scalar + other.scalar,
                vector + other.vector,
                bivector + other.bivector,
                trivector + other.trivector,
                quadvector + other.quadvector
            );
        }

        /**
         * @brief Subtraction operator overload.
         * @param other The multivector to subtract.
         * @return The resulting multivector.
         */
        constexpr Multivector operator-(const Multivector& other) const {
            return Multivector(
                scalar - other.scalar,
                vector - other.vector,
                bivector - other.bivector,
                trivector - other.trivector,
                quadvector - other.quadvector
            );
        }

        /**
         * @brief Unary negation operator overload.
         * @return The resulting multivector.
         */
        constexpr Multivector operator-() const {
            return Multivector(-scalar, -vector, -bivector, -trivector, -quadvector);
        }

        /**
         * @brief Scalar multiplication operator overload.
         * @param s The scalar to multiply by.
         * @return The resulting multivector.
         */
        constexpr Multivector operator*(const Scalar<T>& s) const {
            return Multivector(scalar * s, vector * s, bivector * s, trivector * s, quadvector * s);
        }

        /**
         * @brief Scalar division operator overload.
         * @param s The scalar to divide by.
         * @return The resulting multivector.
         */
        constexpr Multivector operator/(const Scalar<T>& s) const {
            return Multivector(scalar / s, vector / s, bivector / s, trivector / s, quadvector / s);
        }

        /**
         * @brief Adds another multivector in place.
         * @param other The multivector to add.
         * @return A reference to this multivector.
         */
        constexpr Multivector& operator+=(const Multivector& other) {
            *this = *this + other;
            return *this;
        }

        /**
         * @brief Subtracts another multivector in place.
         * @param other The multivector to subtract.
         * @return A reference to this multivector.
         */
        constexpr Multivector& operator-=(const Multivector& other) {
            *this = *this - other;
            return *this;
        }

        /**
         * @brief Scales this multivector in place by a scalar.
         * @param s The scalar to multiply by.
         * @return A reference to this multivector.
         */
        constexpr Multivector& operator*=(const Scalar<T>& s) {
            *this = *this * s;
            return *this;
        }

        /**
         * @brief Divides this multivector in place by a scalar.
         * @param s The scalar to divide by.
         * @return A reference to this multivector.
         */
        constexpr Multivector& operator/=(const Scalar<T>& s) {
            *this = *this / s;
            return *this;
        }

        /**
         * @brief Returns a string representation of the multivector.
         * @return A string representing the multivector.
         */
        std::string to_string() const {
            // Parenthesised so the joining " + " is distinguishable from the
            // "+" inside each component.
            return "(" + scalar.to_string() + ") + (" + vector.to_string()
                 + ") + (" + bivector.to_string() + ") + (" + trivector.to_string()
                 + ") + (" + quadvector.to_string() + ")";
        }
    };

    template<typename T>
    /**
     * @brief Scales a multivector with the raw numeric value on the left.
     * @param value The value to multiply by.
     * @param m The multivector to scale.
     * @return The resulting multivector.
     */
    constexpr Multivector<T> operator*(T value, const Multivector<T>& m) {
        return m * Scalar<T>(value);
    }
} // namespace CliffordCore::PGA
