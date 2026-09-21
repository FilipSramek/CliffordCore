#pragma once

/**
 * @file cliffordcore.hpp
 * @brief Library-wide version information. Not an umbrella -- pick an algebra.
 *
 * This header deliberately includes no algebra. CliffordCore ships one
 * namespace per algebra, and you want exactly one of them:
 *
 * - `#include <cliffordcore/cl3.hpp>` for Cl(3,0), namespace CliffordCore::Cl3
 * - `#include <cliffordcore/pga.hpp>` for Cl(3,0,1), namespace CliffordCore::PGA
 *
 * Cl(2,0) (namespace Cl2) is planned and will slot in the same way. There is no
 * header that pulls in all of them, because mixing two algebras in one
 * expression is a bug, and keeping them in separate namespaces makes it a
 * compile error.
 *
 * @author Filip Sramek
 * @version 0.2.0
 * @date 2026
 * @copyright Copyright (c) 2026 Filip Sramek. All rights reserved.
 * @note Licensing is not settled yet; see README.md.
 */

/**
 * @mainpage CliffordCore
 *
 * CliffordCore is a header-only C++17 library for Clifford (geometric) algebra.
 * Every type is a template over an arithmetic component type, and every
 * operation is constexpr apart from the printing helpers. There is no build
 * system and no dependency -- add `include/` to your include path and go.
 *
 * @section mainpage_algebras One namespace per algebra
 *
 * Each algebra lives in its own namespace under `CliffordCore`, with the **same
 * type names** in each: `Scalar`, `Vector`, `Bivector`, `Multivector`, and so
 * on. Choosing an algebra is therefore one line, and switching is a one-line
 * edit:
 *
 * | Algebra | Header | Namespace | Status |
 * | --- | --- | --- | --- |
 * | Cl(3,0), 3D Euclidean | `<cliffordcore/cl3.hpp>` | `CliffordCore::Cl3` | shipping |
 * | Cl(3,0,1), 3D projective | `<cliffordcore/pga.hpp>` | `CliffordCore::PGA` | shipping |
 * | Cl(2,0), 2D Euclidean | `<cliffordcore/cl2.hpp>` | `CliffordCore::Cl2` | planned |
 *
 * The types are written independently per algebra rather than generated from a
 * signature, because the algebras share a vocabulary rather than an
 * implementation. A grade 1 element is a *direction* in Cl(3,0) but a *plane*
 * in PGA, with a different component count; Cl(2,0)'s bivector is also its
 * pseudoscalar. Names carry across, meanings do not.
 *
 * @section mainpage_usage Usage
 *
 * @code
 * #include <cliffordcore/cl3.hpp>
 * namespace ga = CliffordCore::Cl3;
 *
 * ga::Vector<double> a(1, 2, 3), b(4, 5, 6);
 * auto product = a * b;              // Multivector: dot part + wedge part
 * auto r = ga::rotor_from_axis_angle(ga::Vector<double>(0, 0, 1), 1.5708);
 * auto turned = ga::rotate(a, r);
 * @endcode
 *
 * Compile with `g++ -std=c++17 -Iinclude your_file.cpp`.
 *
 * @section mainpage_basis The Cl(3,0) basis
 *
 * Vectors are `e1 e2 e3` with `e_i^2 = +1`. Bivectors are stored as
 * `(xy, xz, yz) = (e1e2, e1e3, e2e3)` -- note `xz`, not the more common `zx`.
 * The pseudoscalar is `e123 = e1e2e3`. It follows that
 * `e12^2 = e13^2 = e23^2 = e123^2 = -1`, and those minus signs drive every
 * inverse and dual formula in the library. When a result's sign looks wrong,
 * the component ordering is the first thing to check.
 *
 * @section mainpage_pga_basis The Cl(3,0,1) basis
 *
 * Four basis vectors, `e0 e1 e2 e3`, with `e0^2 = 0` and `e_i^2 = +1`
 * otherwise. Every grade is stored lexicographically, matching Cl(3,0)'s
 * choice: `(e01, e02, e03, e12, e13, e23)` and `(e012, e013, e023, e123)`,
 * with the pseudoscalar `e0123`. Grade 1 is a plane, grade 2 a line and
 * grade 3 a point, and the aliases `Plane`, `Line` and `Point` say so.
 *
 * The degenerate `e0` carries an object's offset from the origin, which is
 * what lets a plane, a line or a point sit anywhere in space rather than
 * always passing through zero. It also breaks three Cl(3,0) formulas: the
 * dual becomes a complement rather than a product with the pseudoscalar, a
 * motor's inverse needs a study number, and a norm sees only the components
 * without an `e0` in them. `docs/pga.md` spells all of that out.
 *
 * @section mainpage_more Further reading
 *
 * - `docs/getting-started.md` -- installing and a first program
 * - `docs/conventions.md` -- sign and return-type rules, and the known gaps
 * - `docs/rotations.md` -- rotors in practice
 * - `docs/cheatsheet.md` -- the return-type matrices on one screen
 * - `docs/pga.md` -- the Cl(3,0,1) conventions, which share none of Cl(3,0)'s signs
 * - `examples/` -- six runnable programs, build them with `./build_examples.sh`
 *
 * @section mainpage_includes About the umbrella headers
 *
 * Each algebra's umbrella is a pure aggregate: it declares nothing of its own.
 * Every header it lists is also a valid standalone entry point, so include just
 * the one or two you need if compile time matters.
 *
 * Note that `tests/test_core.cpp` and `tests/test_pga.cpp` deliberately do NOT
 * use an umbrella. Their explicit per-header includes are a live check that
 * each header pulls in its own dependencies; an umbrella include would mask
 * exactly that class of bug, which has broken this codebase twice. Each suite
 * is also its own translation unit, so the two algebras never meet in one.
 */

/**
 * @name Version
 * Header-only libraries get vendored by copy-paste, so the version has to
 * travel with the source rather than living only in a tag.
 * @{
 */
#define CLIFFORDCORE_VERSION_MAJOR 0
#define CLIFFORDCORE_VERSION_MINOR 2
#define CLIFFORDCORE_VERSION_PATCH 0
#define CLIFFORDCORE_VERSION_STRING "0.2.0"
/** @} */
