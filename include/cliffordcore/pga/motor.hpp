#pragma once

/**
 * @file cliffordcore/pga/motor.hpp
 * @brief The even subalgebra: the Motor type, a general rigid motion.
 *
 * Motor<T> holds grades 0, 2 and 4 -- a scalar, all six bivector
 * components, and e0123 -- which is the full even subalgebra of Cl(3,0,1).
 * A motor is a screw motion: a rotation about any line combined with a
 * translation along it, which covers every proper rigid motion of 3D space.
 * Motors compose by multiplication without widening and apply through the
 * sandwich M x ~M to planes, lines, points and whole multivectors.
 *
 * Rotor and Translator widen to Motor implicitly, because nothing is lost.
 * The constructor from Multivector is explicit, because narrowing discards
 * the odd grades.
 *
 * A motor's own product with its reverse is not in general a scalar but a
 * "study number" a + b e0123; inverse.hpp and normalize.hpp account for the
 * b term, and plain division by the scalar part is wrong for a motor that
 * has drifted off unit norm.
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
#include "bivector.hpp"
#include "quadvector.hpp"
#include "rotor.hpp"
#include "translator.hpp"
#include "multivector.hpp"

namespace CliffordCore::PGA
{
    template<typename T>
    /**
     * @brief A class representing a rigid motion: the even subalgebra of Cl(3,0,1).
     * @tparam T The arithmetic component type.
     */
    class Motor
    {
        static_assert(std::is_arithmetic<T>::value, "Motor can only be instantiated with numeric types.");
    public:
        Scalar<T> scalar;          ///< The grade 0 part.
        Bivector<T> bivector;      ///< The grade 2 part: the screw axis, scaled.
        Quadvector<T> quadvector;  ///< The grade 4 part, which couples rotation to translation.

        /**
         * @brief Default constructor initializes all components to zero.
         */
        constexpr Motor() : scalar(), bivector(), quadvector() {}

        /**
         * @brief Constructor initializes the motor with the provided components.
         * @param s The scalar component.
         * @param b The bivector component.
         * @param q The quadvector component.
         */
        constexpr Motor(const Scalar<T>& s, const Bivector<T>& b, const Quadvector<T>& q)
            : scalar(s), bivector(b), quadvector(q) {}

        /**
         * @brief Constructor widens a rotor to a motor. Nothing is lost.
         * @param r The rotor to widen.
         */
        constexpr Motor(const Rotor<T>& r)
            : scalar(r.scalar), bivector(Bivector<T>(0, 0, 0, r.e12, r.e13, r.e23)), quadvector() {}

        /**
         * @brief Constructor widens a translator to a motor. Nothing is lost.
         * @param t The translator to widen.
         */
        constexpr Motor(const Translator<T>& t)
            : scalar(t.scalar), bivector(Bivector<T>(t.e01, t.e02, t.e03, 0, 0, 0)), quadvector() {}

        /**
         * @brief Constructor narrows a multivector to its even part.
         * @param m The multivector to narrow.
         *
         * EXPLICIT on purpose: this discards grades 1 and 3. Spell the
         * narrowing out, or use to_motor(m), which does the same thing by name.
         */
        constexpr explicit Motor(const Multivector<T>& m)
            : scalar(m.scalar), bivector(m.bivector), quadvector(m.quadvector) {}

        /**
         * @brief Addition operator overload.
         * @param other The motor to add.
         * @return The resulting motor.
         */
        constexpr Motor operator+(const Motor& other) const {
            return Motor(scalar + other.scalar, bivector + other.bivector, quadvector + other.quadvector);
        }

        /**
         * @brief Subtraction operator overload.
         * @param other The motor to subtract.
         * @return The resulting motor.
         */
        constexpr Motor operator-(const Motor& other) const {
            return Motor(scalar - other.scalar, bivector - other.bivector, quadvector - other.quadvector);
        }

        /**
         * @brief Unary negation operator overload.
         * @return The resulting motor, which represents the same motion.
         */
        constexpr Motor operator-() const {
            return Motor(-scalar, -bivector, -quadvector);
        }

        /**
         * @brief Scalar multiplication operator overload.
         * @param other The scalar to multiply by.
         * @return The resulting motor.
         */
        constexpr Motor operator*(const Scalar<T>& other) const {
            return Motor(scalar * other, bivector * other, quadvector * other);
        }

        /**
         * @brief Scalar division operator overload.
         * @param other The scalar to divide by.
         * @return The resulting motor.
         */
        constexpr Motor operator/(const Scalar<T>& other) const {
            return Motor(scalar / other, bivector / other, quadvector / other);
        }

        /**
         * @brief Adds another motor in place.
         * @param other The motor to add.
         * @return A reference to this motor.
         */
        constexpr Motor& operator+=(const Motor& other) {
            *this = *this + other;
            return *this;
        }

        /**
         * @brief Subtracts another motor in place.
         * @param other The motor to subtract.
         * @return A reference to this motor.
         */
        constexpr Motor& operator-=(const Motor& other) {
            *this = *this - other;
            return *this;
        }

        /**
         * @brief Scales this motor in place by a scalar.
         * @param other The scalar to multiply by.
         * @return A reference to this motor.
         */
        constexpr Motor& operator*=(const Scalar<T>& other) {
            *this = *this * other;
            return *this;
        }

        /**
         * @brief Divides this motor in place by a scalar.
         * @param other The scalar to divide by.
         * @return A reference to this motor.
         */
        constexpr Motor& operator/=(const Scalar<T>& other) {
            *this = *this / other;
            return *this;
        }

        /**
         * @brief Returns a string representation of the motor.
         * @return A string representing the motor.
         */
        std::string to_string() const {
            return "(" + scalar.to_string() + ") + (" + bivector.to_string()
                 + ") + (" + quadvector.to_string() + ")";
        }
    };

    template<typename T>
    /**
     * @brief Scales a motor with the raw numeric value on the left.
     * @param value The value to multiply by.
     * @param m The motor to scale.
     * @return The resulting motor.
     */
    constexpr Motor<T> operator*(T value, const Motor<T>& m) {
        return m * Scalar<T>(value);
    }
} // namespace CliffordCore::PGA
