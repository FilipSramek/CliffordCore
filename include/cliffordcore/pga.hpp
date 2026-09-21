#pragma once

/**
 * @file cliffordcore/pga.hpp
 * @brief Umbrella header for Cl(3,0,1): every type and operation of the 3D
 *        projective geometric algebra.
 *
 * Including this one header brings in all nine Cl(3,0,1) types and every
 * operation on them, all inside namespace CliffordCore::PGA. It declares
 * nothing itself; each individual header under `cliffordcore/pga/` is also a
 * valid standalone entry point, so include only what you need if compile time
 * matters.
 *
 * The type names are the same as in CliffordCore::Cl3 -- Scalar, Vector,
 * Bivector, Trivector, Multivector, Rotor -- plus Quadvector, Translator and
 * Motor, which Cl(3,0) has no use for. Their meanings differ: here a Vector
 * is a plane, a Bivector a line and a Trivector a point, and the aliases
 * Plane, Line, Point and Twist say so. See docs/pga.md for the conventions.
 *
 * The recommended way to use it is a namespace alias, which makes the choice of
 * algebra one visible line at the top of the file:
 *
 * @code
 * #include <cliffordcore/pga.hpp>
 * namespace ga = CliffordCore::PGA;
 *
 * ga::Point<double> p = ga::point(1.0, 2.0, 3.0);
 * ga::Translator<double> t = ga::translator(0.0, 0.0, 5.0);
 * auto moved = ga::sandwich(p, t);
 * @endcode
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include "../cliffordcore.hpp"

// Types, in grade order. (scalar.hpp alone already pulls in all nine via its
// declare-then-include pattern; the rest are listed for readability.)
#include "pga/scalar.hpp"
#include "pga/vector.hpp"
#include "pga/bivector.hpp"
#include "pga/trivector.hpp"
#include "pga/quadvector.hpp"
#include "pga/multivector.hpp"
#include "pga/rotor.hpp"
#include "pga/translator.hpp"
#include "pga/motor.hpp"

// Products.
#include "pga/operations/dot_product.hpp"
#include "pga/operations/wedge_product.hpp"
#include "pga/operations/geometric_product.hpp"
#include "pga/operations/mixed_products.hpp"

// Sums.
#include "pga/operations/addition.hpp"
#include "pga/operations/subtraction.hpp"

// Magnitudes and involutions.
#include "pga/operations/norm.hpp"
#include "pga/operations/normalize.hpp"
#include "pga/operations/reverse.hpp"
#include "pga/operations/involutions.hpp"
#include "pga/operations/inverse.hpp"

// Structure.
#include "pga/operations/dual.hpp"
#include "pga/operations/regressive_product.hpp"
#include "pga/operations/grade.hpp"
#include "pga/operations/contraction.hpp"
#include "pga/operations/comparison.hpp"
#include "pga/operations/stream.hpp"
#include "pga/operations/geometry.hpp"
#include "pga/operations/primitives.hpp"

// Motions.
#include "pga/operations/exp.hpp"
#include "pga/operations/log.hpp"
#include "pga/operations/sandwich.hpp"
#include "pga/operations/motor_construction.hpp"
