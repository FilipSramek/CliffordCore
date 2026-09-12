#pragma once

/**
 * @file cliffordcore/cl3.hpp
 * @brief Umbrella header for Cl(3,0): every type and operation of the 3D
 *        Euclidean geometric algebra.
 *
 * Including this one header brings in all six Cl(3,0) types and every
 * operation on them, all inside namespace CliffordCore::Cl3. It declares
 * nothing itself; each individual header under `cliffordcore/cl3/` is also a
 * valid standalone entry point, so include only what you need if compile time
 * matters.
 *
 * The recommended way to use it is a namespace alias, which makes the choice of
 * algebra one visible line at the top of the file:
 *
 * @code
 * #include <cliffordcore/cl3.hpp>
 * namespace ga = CliffordCore::Cl3;
 *
 * ga::Vector<double> a(1, 2, 3), b(4, 5, 6);
 * auto product = a * b;
 * @endcode
 *
 * @author Filip Sramek
 * @version 0.1.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

#include "../cliffordcore.hpp"

// Types, in grade order. (scalar.hpp alone already pulls in all six via its
// declare-then-include pattern; the rest are listed for readability.)
#include "cl3/scalar.hpp"
#include "cl3/vector.hpp"
#include "cl3/bivector.hpp"
#include "cl3/trivector.hpp"
#include "cl3/multivector.hpp"
#include "cl3/rotor.hpp"

// Products.
#include "cl3/operations/dot_product.hpp"
#include "cl3/operations/wedge_product.hpp"
#include "cl3/operations/contraction.hpp"
#include "cl3/operations/geometric_product.hpp"
#include "cl3/operations/mixed_products.hpp"

// Sums.
#include "cl3/operations/addition.hpp"
#include "cl3/operations/subtraction.hpp"

// Magnitudes and involutions.
#include "cl3/operations/norm.hpp"
#include "cl3/operations/normalize.hpp"
#include "cl3/operations/reverse.hpp"
#include "cl3/operations/involutions.hpp"
#include "cl3/operations/inverse.hpp"

// Structure.
#include "cl3/operations/dual.hpp"
#include "cl3/operations/grade.hpp"
#include "cl3/operations/comparison.hpp"
#include "cl3/operations/stream.hpp"
#include "cl3/operations/geometry.hpp"

// Rotations.
#include "cl3/operations/exp.hpp"
#include "cl3/operations/log.hpp"
#include "cl3/operations/sandwich.hpp"
#include "cl3/operations/rotor_construction.hpp"
