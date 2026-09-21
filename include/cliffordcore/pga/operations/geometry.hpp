#pragma once

/**
 * @file cliffordcore/pga/operations/geometry.hpp
 * @brief The classical answers: reflect, project, distance, angle.
 *
 * The questions a 3D geometry library is actually asked. In PGA they are
 * one-liners over the products, and -- unlike the vector-algebra versions
 * -- the same formula covers a line that misses the origin and a plane
 * anywhere in space, because the offset rides in the e0 components.
 *
 * Projection is (x _| y) * inverse(y), the same shape as Cl3's. Reflection
 * in a plane is p x ~p with a grade-dependent sign: planes and lines flip,
 * points do not, which keeps a reflected point's weight positive.
 *
 * The distance and angle functions divide by norms, so they return a raw
 * arithmetic value rather than a Scalar, and they use std::acos and
 * std::sqrt, which makes them inline rather than constexpr -- the same
 * exception approx_equal makes.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include <cmath>

#include "../scalar.hpp"
#include "../vector.hpp"
#include "../bivector.hpp"
#include "../trivector.hpp"
#include "../multivector.hpp"
#include "contraction.hpp"
#include "dot_product.hpp"
#include "geometric_product.hpp"
#include "grade.hpp"
#include "inverse.hpp"
#include "mixed_products.hpp"
#include "norm.hpp"
#include "regressive_product.hpp"
#include "wedge_product.hpp"

namespace CliffordCore::PGA
{
    // -----------------------------------------------------------------------
    // Reflection in a plane
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Reflects a plane in another plane.
     * @param x The plane to reflect.
     * @param mirror The plane to reflect in. Need not be normalized.
     * @return The reflected plane, -p x p^-1.
     */
    constexpr Vector<T> reflect(const Vector<T>& x, const Vector<T>& mirror) {
        return -grade1(mirror * x * inverse(mirror));
    }

    template<typename T>
    /**
     * @brief Reflects a line in a plane.
     * @param x The line to reflect.
     * @param mirror The plane to reflect in. Need not be normalized.
     * @return The reflected line, -p x p^-1. A line lying in the mirror is
     *         unchanged up to orientation; one perpendicular to it is reversed.
     */
    constexpr Bivector<T> reflect(const Bivector<T>& x, const Vector<T>& mirror) {
        return -grade2(mirror * x * inverse(mirror));
    }

    template<typename T>
    /**
     * @brief Reflects a point in a plane.
     * @param x The point to reflect.
     * @param mirror The plane to reflect in. Need not be normalized.
     * @return The reflected point, p x p^-1. Points take the opposite sign to
     *         planes and lines here, which is what keeps a reflected point's
     *         weight positive rather than flipping it.
     */
    constexpr Trivector<T> reflect(const Trivector<T>& x, const Vector<T>& mirror) {
        return grade3(mirror * x * inverse(mirror));
    }

    // -----------------------------------------------------------------------
    // Orthogonal projection
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief Projects a point onto a plane.
     * @param x The point to project.
     * @param onto The plane to project onto. Need not be normalized.
     * @return The closest point of the plane to x: (onto _| x) * onto,
     *         normalized by the plane's squared norm.
     */
    constexpr Trivector<T> project(const Trivector<T>& x, const Vector<T>& onto) {
        const T n = squared_norm(onto).value;
        if (n == T(0)) {
            return Trivector<T>();
        }
        return grade3(left_contraction(onto, x) * onto) / Scalar<T>(n);
    }

    template<typename T>
    /**
     * @brief Projects a point onto a line.
     * @param x The point to project.
     * @param onto The line to project onto. Need not be normalized.
     * @return The closest point of the line to x. An ideal line has no
     *         projection; zero is returned for one.
     */
    constexpr Trivector<T> project(const Trivector<T>& x, const Bivector<T>& onto) {
        const T n = squared_norm(onto).value;
        if (n == T(0)) {
            return Trivector<T>();
        }
        // inverse(onto) is -onto / |onto|^2 for a line, which supplies the sign.
        return grade3(left_contraction(onto, x) * inverse(onto));
    }

    template<typename T>
    /**
     * @brief Projects a line onto a plane.
     * @param x The line to project.
     * @param onto The plane to project onto. Need not be normalized.
     * @return The line in the plane closest to x. An ideal plane gives zero.
     */
    constexpr Bivector<T> project(const Bivector<T>& x, const Vector<T>& onto) {
        const T n = squared_norm(onto).value;
        if (n == T(0)) {
            return Bivector<T>();
        }
        return grade2(left_contraction(onto, x) * onto) / Scalar<T>(n);
    }

    // -----------------------------------------------------------------------
    // Distances
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief The Euclidean distance between two points.
     * @param a The first point.
     * @param b The second point.
     * @return norm(join(a, b)) divided by the two weights. Zero is returned if
     *         either point is ideal, since a point at infinity has no distance
     *         to anything.
     */
    inline T distance(const Trivector<T>& a, const Trivector<T>& b) {
        const T wa = a.e123;
        const T wb = b.e123;
        if (wa == T(0) || wb == T(0)) {
            return T(0);
        }
        return norm(join(a, b)).value / std::abs(wa * wb);
    }

    template<typename T>
    /**
     * @brief The signed distance from a point to a plane.
     * @param p The point.
     * @param pl The plane.
     * @return (pl ^ p) divided by the plane's norm and the point's weight. The
     *         sign says which side of the plane the point is on; it is positive
     *         on the side the normal points away from. Zero if the point is
     *         ideal or the plane has no normal.
     */
    inline T distance(const Trivector<T>& p, const Vector<T>& pl) {
        const T w = p.e123;
        const T n = norm(pl).value;
        if (w == T(0) || n == T(0)) {
            return T(0);
        }
        return (pl ^ p).e0123 / (n * w);
    }

    template<typename T>
    /**
     * @brief The distance from a point to a line.
     * @param p The point.
     * @param l The line.
     * @return norm(join(p, l)) divided by the point's weight and the line's
     *         norm: the join is the plane through both, and its Euclidean norm
     *         is the perpendicular distance. Zero if the point is ideal or the
     *         line is.
     */
    inline T distance(const Trivector<T>& p, const Bivector<T>& l) {
        const T w = p.e123;
        const T n = norm(l).value;
        if (w == T(0) || n == T(0)) {
            return T(0);
        }
        return norm(join(p, l)).value / std::abs(w * n);
    }

    // -----------------------------------------------------------------------
    // Angles
    // -----------------------------------------------------------------------

    template<typename T>
    /**
     * @brief The angle between two planes, in radians.
     * @param a The first plane.
     * @param b The second plane.
     * @return acos of their normalized dot product, in [0, pi]. Zero if either
     *         plane is ideal. The planes' positions do not matter -- only the
     *         normals, which is what the degenerate metric arranges.
     */
    inline T angle(const Vector<T>& a, const Vector<T>& b) {
        const T na = norm(a).value;
        const T nb = norm(b).value;
        if (na == T(0) || nb == T(0)) {
            return T(0);
        }
        const T c = (a | b).value / (na * nb);
        return std::acos(c < T(-1) ? T(-1) : (c > T(1) ? T(1) : c));
    }

    template<typename T>
    /**
     * @brief The angle between two lines, in radians.
     * @param a The first line.
     * @param b The second line.
     * @return acos of minus their normalized scalar product, in [0, pi]. The
     *         minus sign is there because a unit line squares to -1. Zero if
     *         either line is ideal; skew lines still give the angle between
     *         their directions.
     */
    inline T angle(const Bivector<T>& a, const Bivector<T>& b) {
        const T na = norm(a).value;
        const T nb = norm(b).value;
        if (na == T(0) || nb == T(0)) {
            return T(0);
        }
        const T c = -grade0(a * b).value / (na * nb);
        return std::acos(c < T(-1) ? T(-1) : (c > T(1) ? T(1) : c));
    }
} // namespace CliffordCore::PGA
