#pragma once

/**
 * @file cliffordcore/cl2.hpp
 * @brief Umbrella header for Cl(2,0): every type and operation of the 2D
 *        Euclidean geometric algebra.
 *
 * Including this one header brings in all five Cl(2,0) types and every
 * operation on them, all inside namespace CliffordCore::Cl2. It declares
 * nothing itself; each individual header under `cliffordcore/cl2/` is also a
 * valid standalone entry point, so include only what you need if compile time
 * matters.
 *
 * The recommended way to use it is a namespace alias, which makes the choice of
 * algebra one visible line at the top of the file:
 *
 * @code
 * #include <cliffordcore/cl2.hpp>
 * namespace ga = CliffordCore::Cl2;
 *
 * ga::Vector<double> a(1, 2), b(3, 4);
 * auto product = a * b;
 * @endcode
 *
 * @author Filip Sramek
 * @version 0.3.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. Licensed under the Apache License, Version 2.0.
 * @note SPDX-License-Identifier: Apache-2.0 -- see the LICENSE file for the full text.
 */

#include "../cliffordcore.hpp"

// Types, in grade order. (scalar.hpp alone already pulls in all five via its
// declare-then-include pattern; the rest are listed for readability.)
#include "cl2/scalar.hpp"
#include "cl2/vector.hpp"
#include "cl2/bivector.hpp"
#include "cl2/multivector.hpp"
#include "cl2/rotor.hpp"
