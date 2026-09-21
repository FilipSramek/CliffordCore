#pragma once

/**
 * @file cliffordcore/pga/operations/motor_construction.hpp
 * @brief Ways to build a rotor, a translator or a motor without assembling components by hand.
 *
 * Orientation convention, the same as Cl3's: exp(theta B) turns by 2 theta
 * in the plane B, and sandwich() applies M x ~M. Working that through, a
 * rotation by `angle` about a unit line L is exp(-(angle/2) L), which is
 * the ordinary right-hand rule -- +90 degrees about +z takes point(1,0,0)
 * to point(0,1,0), and plane x=0 to plane y=0, exactly as Cl3's
 * rotor_from_axis_angle does for vectors.
 *
 * What PGA adds is that the line need not pass through the origin, and
 * that a translation along the axis can ride along: screw(L, angle, d) is
 * the general rigid motion, and every motor is one.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include <cmath>
#include <limits>

#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../quadvector.hpp"
#include "../rotor.hpp"
#include "../translator.hpp"
#include "../motor.hpp"
#include "exp.hpp"
#include "log.hpp"
#include "geometric_product.hpp"
#include "mixed_products.hpp"
#include "norm.hpp"
#include "normalize.hpp"
#include "primitives.hpp"
#include "reverse.hpp"

namespace CliffordCore::PGA
{
    // -----------------------------------------------------------------------
    // Identities
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief The rotor that represents no rotation.
     * @return The identity rotor, scalar 1 with a zero bivector.
     */
    constexpr Rotor<T> identity_rotor() {
        return Rotor<T>(Scalar<T>(1), T(0), T(0), T(0));
    }

    template<typename T>
    /**
     * @brief The translator that represents no translation.
     * @return The identity translator, scalar 1 with a zero ideal part.
     */
    constexpr Translator<T> identity_translator() {
        return Translator<T>(Scalar<T>(1), T(0), T(0), T(0));
    }

    template<typename T>
    /**
     * @brief The motor that represents no motion.
     * @return The identity motor, scalar 1 with everything else zero.
     */
    constexpr Motor<T> identity_motor() {
        return Motor<T>(Scalar<T>(1), Bivector<T>(), Quadvector<T>());
    }

    // -----------------------------------------------------------------------
    // Translators
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Builds the translator that displaces by (dx, dy, dz).
     * @param dx The x displacement.
     * @param dy The y displacement.
     * @param dz The z displacement.
     * @return 1 - (dx e01 + dy e02 + dz e03) / 2. The half is what the
     *         two-sided sandwich needs; the minus sign is fixed by the point
     *         convention in primitives.hpp.
     */
    constexpr Translator<T> translator(T dx, T dy, T dz) {
        return Translator<T>(Scalar<T>(1), -dx / T(2), -dy / T(2), -dz / T(2));
    }

    template<typename T>
    /**
     * @brief Builds the translator taking one point to another.
     * @param from The starting point.
     * @param to The target point.
     * @return A normalised translator t with sandwich(from, t) equal to to.
     *         Points of zero weight have no position, so the identity is
     *         returned if either is ideal.
     */
    constexpr Translator<T> translator_between(const Trivector<T>& from, const Trivector<T>& to) {
        if (from.e123 == T(0) || to.e123 == T(0)) {
            return identity_translator<T>();
        }
        const std::array<T, 3> a = position(from);
        const std::array<T, 3> b = position(to);
        return translator<T>(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
    }

    // -----------------------------------------------------------------------
    // Rotations about a line
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Builds a rotor from an axis direction and an angle.
     * @param dx The x component of the axis direction.
     * @param dy The y component of the axis direction.
     * @param dz The z component of the axis direction.
     * @param angle The rotation angle in radians, right-handed about the axis.
     * @return The rotor turning by that angle about the line through the
     *         ORIGIN in that direction, or the identity if the direction is
     *         zero. For an axis that misses the origin, use
     *         motor_from_line_angle with a line built by line_through_points.
     */
    constexpr Rotor<T> rotor_from_axis_angle(T dx, T dy, T dz, T angle) {
        const T length = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (length == T(0)) {
            return identity_rotor<T>();
        }
        const T half = angle / T(2);
        const T s = std::sin(half) / length;
        // exp(-(angle/2) L) with L = line_through_origin(dx, dy, dz),
        // whose components are (e12, e13, e23) = (dz, -dy, dx).
        return Rotor<T>(Scalar<T>(std::cos(half)), -s * dz, s * dy, -s * dx);
    }

    template<typename T>
    /**
     * @brief Builds a motor that rotates about an arbitrary line.
     * @param axis The line to rotate about. It is normalized internally, so any
     *        non-ideal line works; the line need not pass through the origin.
     * @param angle The rotation angle in radians, right-handed about the line's direction.
     * @return exp(-(angle/2) * normalize(axis)), or the identity motor if the
     *         line is ideal.
     */
    constexpr Motor<T> motor_from_line_angle(const Bivector<T>& axis, T angle) {
        if (norm(axis).value == T(0)) {
            return identity_motor<T>();
        }
        return exp(normalize(axis) * Scalar<T>(-angle / T(2)));
    }

    template<typename T>
    /**
     * @brief Builds a screw motion: a rotation about a line plus a translation along it.
     * @param axis The screw axis. Normalized internally; it need not pass through the origin.
     * @param angle The rotation angle in radians, right-handed about the line's direction.
     * @param distance The distance to travel along the axis while turning.
     * @return The motor performing both at once. Chasles' theorem says every
     *         rigid motion is one of these. An ideal axis gives the identity.
     */
    constexpr Motor<T> screw(const Bivector<T>& axis, T angle, T distance) {
        if (norm(axis).value == T(0)) {
            return identity_motor<T>();
        }
        const Bivector<T> l = normalize(axis);
        // l * e0123 is the ideal bivector that translates along l:
        // (e12, e13, e23) map onto (-e03, +e02, -e01).
        const Bivector<T> along(-l.e23, l.e13, -l.e12, T(0), T(0), T(0));
        return exp(l * Scalar<T>(-angle / T(2)) + along * Scalar<T>(distance / T(2)));
    }

    // -----------------------------------------------------------------------
    // Motors from the geometry they should produce
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Builds the motor that takes one plane onto another.
     * @param from The starting plane. Normalized internally.
     * @param to The target plane. Normalized internally.
     * @return A unit motor m with sandwich(from, m) equal to to: the rotation
     *         about their common line by the angle between them, or a
     *         translation if they are parallel.
     *
     * Built as normalize(1 + to*from), the half-angle construction -- the same
     * shape as Cl3's rotor_between. When the planes are opposite the motion is
     * ambiguous; a half turn about a definitely non-degenerate line is chosen.
     */
    constexpr Motor<T> motor_between(const Vector<T>& from, const Vector<T>& to) {
        if (norm(from).value == T(0) || norm(to).value == T(0)) {
            return identity_motor<T>();
        }

        const Vector<T> a = normalize(from);
        const Vector<T> b = normalize(to);
        const T alignment = (a | b).value;

        // Opposite normals: every line in `a` gives a valid half turn, so pick
        // one that is definitely not degenerate.
        if (alignment < T(-1) + std::numeric_limits<T>::epsilon() * T(8)) {
            const Vector<T> reference = (std::abs(a.e1) < T(0.9))
                ? plane<T>(T(1), T(0), T(0), T(0))
                : plane<T>(T(0), T(1), T(0), T(0));
            return Motor<T>(Scalar<T>(0), normalize(a ^ reference), Quadvector<T>());
        }

        const Motor<T> product = motor_product(b, a);
        return normalize(Motor<T>(product.scalar + Scalar<T>(1), product.bivector, product.quadvector));
    }

    template<typename T>
    /**
     * @brief Combines a rotation about the origin with a translation.
     * @param r The rotation, applied first.
     * @param t The translation, applied second.
     * @return The motor t * r, which rotates and then translates.
     */
    constexpr Motor<T> motor_from_rotor_translator(const Rotor<T>& r, const Translator<T>& t) {
        return t * r;
    }

    // -----------------------------------------------------------------------
    // Interpolation
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Interpolates between two motors along the shortest screw.
     * @param from The motor at t = 0.
     * @param to The motor at t = 1.
     * @param t The interpolation parameter, normally in [0, 1].
     * @return The interpolated motor, a * exp(t * log(reverse(a) * b)).
     *
     * A motor and its negation describe the same motion, so if the two point
     * away from each other the second is negated first and the result takes
     * the short way round. Interpolating this way moves along a single screw:
     * the rotation and the translation stay coupled, which is what makes it
     * the right thing for animating a rigid body.
     */
    constexpr Motor<T> slerp(const Motor<T>& from, const Motor<T>& to, T t) {
        const Motor<T> a = normalize(from);
        Motor<T> b = normalize(to);

        const T alignment = a.scalar.value * b.scalar.value
                          + a.bivector.e12 * b.bivector.e12
                          + a.bivector.e13 * b.bivector.e13
                          + a.bivector.e23 * b.bivector.e23;
        if (alignment < T(0)) {
            b = -b;
        }

        // reverse(a) is the inverse of a unit motor, so this is the relative
        // motion from a to b, scaled by t and reapplied.
        return a * exp(log(reverse(a) * b) * Scalar<T>(t));
    }

    template<typename T>
    /**
     * @brief Interpolates between two rotors along the shortest arc.
     * @param from The rotor at t = 0.
     * @param to The rotor at t = 1.
     * @param t The interpolation parameter, normally in [0, 1].
     * @return The interpolated rotor. The same construction as the motor
     *         overload, narrowed back to a Rotor since no translation can
     *         appear between two rotations about the origin.
     */
    constexpr Rotor<T> slerp(const Rotor<T>& from, const Rotor<T>& to, T t) {
        return Rotor<T>(slerp(Motor<T>(from), Motor<T>(to), t));
    }

    template<typename T>
    /**
     * @brief Interpolates between two translators.
     * @param from The translator at t = 0.
     * @param to The translator at t = 1.
     * @param t The interpolation parameter, normally in [0, 1].
     * @return The interpolated translator, which moves along a straight line
     *         at constant speed.
     */
    constexpr Translator<T> slerp(const Translator<T>& from, const Translator<T>& to, T t) {
        return Translator<T>(slerp(Motor<T>(from), Motor<T>(to), t));
    }
} // namespace CliffordCore::PGA
