#pragma once

/**
 * @file cliffordcore/pga/operations/primitives.hpp
 * @brief Building planes, lines and points from ordinary numbers, and reading them back.
 *
 * The component signs of a PGA point are not memorable, and writing them
 * out by hand is the easiest way to get a sign wrong. These factories are
 * the intended way in: point(x, y, z) and plane(a, b, c, d) hide the
 * convention, and position(), normal() and direction() read it back. The
 * conventions themselves are in docs/pga.md.
 *
 * The two ways to build a line are here too: through two points (join) and
 * as the intersection of two planes (meet). Both are named after what they
 * produce rather than after the product they use.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <array>

#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../quadvector.hpp"
#include "regressive_product.hpp"
#include "wedge_product.hpp"

namespace CliffordCore::PGA
{
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Builds the plane a x + b y + c z + d = 0.
     * @param a The x coefficient of the normal.
     * @param b The y coefficient of the normal.
     * @param c The z coefficient of the normal.
     * @param d The offset: the plane's distance from the origin, negated and
     *        scaled by the length of the normal.
     * @return The vector d e0 + a e1 + b e2 + c e3.
     */
    constexpr Vector<T> plane(T a, T b, T c, T d) {
        return Vector<T>(d, a, b, c);
    }

    template<typename T>
    /**
     * @brief Builds the Euclidean point (x, y, z), with weight one.
     * @param x The x coordinate.
     * @param y The y coordinate.
     * @param z The z coordinate.
     * @return The trivector e123 - x e023 + y e013 - z e012. This is the
     *         triple meet of the planes X = x, Y = y and Z = z, so a point
     *         built this way lies on a plane exactly when plane ^ point is zero.
     */
    constexpr Trivector<T> point(T x, T y, T z) {
        return Trivector<T>(-z, y, -x, T(1));
    }

    template<typename T>
    /**
     * @brief Builds the ideal point in the direction (dx, dy, dz): a point at infinity.
     * @param dx The x component of the direction.
     * @param dy The y component of the direction.
     * @param dz The z component of the direction.
     * @return The trivector with the same ideal components as point(dx, dy, dz)
     *         but zero weight. Adding an ideal point to a Euclidean one
     *         translates it; translating an ideal point leaves it unchanged.
     */
    constexpr Trivector<T> ideal_point(T dx, T dy, T dz) {
        return Trivector<T>(-dz, dy, -dx, T(0));
    }

    template<typename T>
    /**
     * @brief The origin.
     * @return The trivector e123, which is point(0, 0, 0).
     */
    constexpr Trivector<T> origin() {
        return Trivector<T>(T(0), T(0), T(0), T(1));
    }

    template<typename T>
    /**
     * @brief The plane at infinity.
     * @return The vector e0, the one plane every ideal point lies on.
     */
    constexpr Vector<T> ideal_plane() {
        return Vector<T>(T(1), T(0), T(0), T(0));
    }

    template<typename T>
    /**
     * @brief Builds the line through the origin in the direction (dx, dy, dz).
     * @param dx The x component of the direction.
     * @param dy The y component of the direction.
     * @param dz The z component of the direction.
     * @return The bivector dz e12 - dy e13 + dx e23, which has no ideal part.
     *         The x, y and z axes are e23, -e13 and e12.
     */
    constexpr Bivector<T> line_through_origin(T dx, T dy, T dz) {
        return Bivector<T>(T(0), T(0), T(0), dz, -dy, dx);
    }

    template<typename T>
    /**
     * @brief Builds the line through two points.
     * @param a The first point.
     * @param b The second point.
     * @return join(a, b). For unit-weight points its Euclidean part is the
     *         direction from a to b and its norm their distance; it is zero if
     *         the points coincide.
     */
    constexpr Bivector<T> line_through_points(const Trivector<T>& a, const Trivector<T>& b) {
        return join(a, b);
    }

    template<typename T>
    /**
     * @brief Builds the line where two planes intersect.
     * @param p The first plane.
     * @param q The second plane.
     * @return p ^ q. Parallel planes give an ideal line; coincident planes give zero.
     */
    constexpr Bivector<T> line_from_planes(const Vector<T>& p, const Vector<T>& q) {
        return p ^ q;
    }

    template<typename T>
    /**
     * @brief Builds the plane through three points.
     * @param a The first point.
     * @param b The second point.
     * @param c The third point.
     * @return join(join(a, b), c). Zero if the three points are collinear.
     */
    constexpr Vector<T> plane_through_points(const Trivector<T>& a, const Trivector<T>& b, const Trivector<T>& c) {
        return join(join(a, b), c);
    }

    // -----------------------------------------------------------------------
    // Reading the geometry back out
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Recovers the Euclidean coordinates of a point.
     * @param p The trivector.
     * @return {x, y, z}, dividing out the weight. An ideal point has no
     *         position; this returns zeroes for one rather than dividing by zero.
     */
    constexpr std::array<T, 3> position(const Trivector<T>& p) {
        if (p.e123 == T(0)) {
            return std::array<T, 3>{T(0), T(0), T(0)};
        }
        return std::array<T, 3>{-p.e023 / p.e123, p.e013 / p.e123, -p.e012 / p.e123};
    }

    template<typename T>
    /**
     * @brief Recovers the direction an ideal point represents, or the unweighted offset of a Euclidean one.
     * @param p The trivector.
     * @return {dx, dy, dz} straight from the ideal components, without dividing
     *         by the weight.
     */
    constexpr std::array<T, 3> direction(const Trivector<T>& p) {
        return std::array<T, 3>{-p.e023, p.e013, -p.e012};
    }

    template<typename T>
    /**
     * @brief Recovers the direction of a line.
     * @param l The bivector.
     * @return {dx, dy, dz} from the Euclidean components: {e23, -e13, e12}.
     *         Zero for an ideal line.
     */
    constexpr std::array<T, 3> direction(const Bivector<T>& l) {
        return std::array<T, 3>{l.e23, -l.e13, l.e12};
    }

    template<typename T>
    /**
     * @brief Recovers the moment of a line.
     * @param l The bivector.
     * @return {e01, e02, e03}, the ideal half of the line.
     */
    constexpr std::array<T, 3> moment(const Bivector<T>& l) {
        return std::array<T, 3>{l.e01, l.e02, l.e03};
    }

    template<typename T>
    /**
     * @brief Recovers the normal of a plane.
     * @param p The vector.
     * @return {a, b, c} = {e1, e2, e3}. Zero for the ideal plane.
     */
    constexpr std::array<T, 3> normal(const Vector<T>& p) {
        return std::array<T, 3>{p.e1, p.e2, p.e3};
    }

    template<typename T>
    /**
     * @brief Recovers the offset of a plane.
     * @param p The vector.
     * @return The e0 component, which is the d of a x + b y + c z + d = 0.
     */
    constexpr Scalar<T> offset(const Vector<T>& p) {
        return Scalar<T>(p.e0);
    }

    // -----------------------------------------------------------------------
    // Classification
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Tests whether a plane is the ideal plane at infinity.
     * @param p The vector.
     * @return True if its normal is zero.
     */
    constexpr bool is_ideal(const Vector<T>& p) {
        return p.e1 == T(0) && p.e2 == T(0) && p.e3 == T(0);
    }

    template<typename T>
    /**
     * @brief Tests whether a line lies at infinity.
     * @param l The bivector.
     * @return True if its direction is zero.
     */
    constexpr bool is_ideal(const Bivector<T>& l) {
        return l.e12 == T(0) && l.e13 == T(0) && l.e23 == T(0);
    }

    template<typename T>
    /**
     * @brief Tests whether a point lies at infinity: a direction rather than a position.
     * @param p The trivector.
     * @return True if its weight is zero.
     */
    constexpr bool is_ideal(const Trivector<T>& p) {
        return p.e123 == T(0);
    }
} // namespace CliffordCore::PGA
