#pragma once

/**
 * @file cliffordcore.hpp
 * @brief Umbrella header: includes every CliffordCore type and operation.
 */

/**
 * @mainpage CliffordCore
 *
 * CliffordCore is a header-only C++17 library for 3D Clifford (geometric)
 * algebra, Cl(3,0). Everything lives in namespace CliffordCore, every type is a
 * template over an arithmetic component type, and every operation is constexpr.
 * There is no build system and no dependency -- add `include/` to your include
 * path and go.
 *
 * @section mainpage_usage Usage
 *
 * @code
 * #include <cliffordcore.hpp>
 *
 * CliffordCore::Vector3<double> a(1, 2, 3), b(4, 5, 6);
 * auto product = a * b;              // Multivector3: dot part + wedge part
 * auto r = CliffordCore::rotor_from_axis_angle(CliffordCore::Vector3<double>(0, 0, 1), 1.5708);
 * auto turned = CliffordCore::rotate(a, r);
 * @endcode
 *
 * Compile with `g++ -std=c++17 -Iinclude your_file.cpp`.
 *
 * @section mainpage_basis The basis
 *
 * Vectors are `e1 e2 e3` with `e_i^2 = +1`. Bivectors are stored as
 * `(xy, xz, yz) = (e1e2, e1e3, e2e3)` -- note `xz`, not the more common `zx`.
 * The pseudoscalar is `e123 = e1e2e3`. It follows that
 * `e12^2 = e13^2 = e23^2 = e123^2 = -1`, and those minus signs drive every
 * inverse and dual formula in the library. When a result's sign looks wrong,
 * the component ordering is the first thing to check.
 *
 * @section mainpage_more Further reading
 *
 * - `docs/getting-started.md` -- installing and a first program
 * - `docs/conventions.md` -- sign and return-type rules, and the known gaps
 * - `docs/rotations.md` -- rotors in practice
 * - `docs/cheatsheet.md` -- the return-type matrices on one screen
 * - `examples/` -- five runnable programs, build them with `./build_examples.sh`
 *
 * @section mainpage_includes About this header
 *
 * This header is a pure aggregate: it declares nothing of its own. Every header
 * it lists is also a valid standalone entry point, so include just the one or
 * two you need if compile time matters.
 *
 * Note that `tests/test_core.cpp` deliberately does NOT use this header. Its
 * explicit per-header includes are a live check that each header pulls in its
 * own dependencies; an umbrella include would mask exactly that class of bug,
 * which has broken this codebase twice.
 */

// Types, in grade order. (scalar.hpp alone already pulls in all six via its
// declare-then-include pattern; the rest are listed for readability.)
#include "scalar.hpp"
#include "vector3.hpp"
#include "bivector3.hpp"
#include "trivector3.hpp"
#include "multivector3.hpp"
#include "rotor3.hpp"

// Products.
#include "operations/dot_product.hpp"
#include "operations/wedge_product.hpp"
#include "operations/contraction.hpp"
#include "operations/geometric_product.hpp"
#include "operations/mixed_products.hpp"

// Sums.
#include "operations/addition.hpp"
#include "operations/subtraction.hpp"

// Magnitudes and involutions.
#include "operations/norm.hpp"
#include "operations/normalize.hpp"
#include "operations/reverse.hpp"
#include "operations/involutions.hpp"
#include "operations/inverse.hpp"

// Structure.
#include "operations/dual.hpp"
#include "operations/grade.hpp"

// Rotations.
#include "operations/exp.hpp"
#include "operations/log.hpp"
#include "operations/sandwich.hpp"
#include "operations/rotor_construction.hpp"
