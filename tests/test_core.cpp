// CliffordCore test suite.
//
// Build and run:
//   g++ -std=c++17 -Iinclude tests/test_core.cpp -o build/test_core && ./build/test_core.exe
//   ...or just ./build.sh
//
// Structure: a small assertion harness, then one section per type and per
// operation. Every check records a pass/fail and the run continues, so a single
// build reports every problem rather than stopping at the first.

#include <cmath>
#include <iostream>
#include <string>

#include "../include/scalar.hpp"
#include "../include/vector3.hpp"
#include "../include/bivector3.hpp"
#include "../include/trivector3.hpp"
#include "../include/multivector3.hpp"
#include "../include/rotor3.hpp"
#include "../include/operations/dot_product.hpp"
#include "../include/operations/wedge_product.hpp"
#include "../include/operations/geometric_product.hpp"
#include "../include/operations/norm.hpp"
#include "../include/operations/inverse.hpp"
#include "../include/operations/reverse.hpp"
#include "../include/operations/dual.hpp"
#include "../include/operations/exp.hpp"
#include "../include/operations/log.hpp"
#include "../include/operations/normalize.hpp"
#include "../include/operations/sandwich.hpp"
#include "../include/operations/addition.hpp"
#include "../include/operations/subtraction.hpp"
#include "../include/operations/mixed_products.hpp"
#include "../include/operations/rotor_construction.hpp"
#include "../include/operations/grade.hpp"

using CliffordCore::Bivector3;
using CliffordCore::Multivector3;
using CliffordCore::Rotor3;
using CliffordCore::Scalar;
using CliffordCore::Trivector3;
using CliffordCore::Vector3;

namespace {

// ---------------------------------------------------------------------------
// Assertion harness
// ---------------------------------------------------------------------------

int g_checks = 0;
int g_failures = 0;
const char* g_section = "";

void section(const char* name)
{
    g_section = name;
}

void fail(const std::string& what, const std::string& detail)
{
    ++g_failures;
    std::cerr << "FAIL [" << g_section << "] " << what;
    if (!detail.empty()) {
        std::cerr << " -- " << detail;
    }
    std::cerr << "\n";
}

// Exact boolean check, for things like "this value is identically zero".
void check(bool condition, const std::string& what)
{
    ++g_checks;
    if (!condition) {
        fail(what, "");
    }
}

// Floating-point comparison. Never compare doubles with ==; values that are
// mathematically equal can differ in the last bits after arithmetic.
void check_close(double actual, double expected, const std::string& what, double tolerance = 1e-12)
{
    ++g_checks;
    if (std::fabs(actual - expected) > tolerance) {
        fail(what, "expected " + std::to_string(expected) + ", got " + std::to_string(actual));
    }
}

// Whole-object helpers keep per-component noise out of the test bodies.
void check_scalar(const Scalar<double>& s, double value, const std::string& what,
                  double tolerance = 1e-12)
{
    check_close(s.value, value, what + ".value", tolerance);
}

void check_vector(const Vector3<double>& v, double x, double y, double z, const std::string& what,
                  double tolerance = 1e-12)
{
    check_close(v.x, x, what + ".x", tolerance);
    check_close(v.y, y, what + ".y", tolerance);
    check_close(v.z, z, what + ".z", tolerance);
}

void check_bivector(const Bivector3<double>& b, double xy, double xz, double yz, const std::string& what,
                    double tolerance = 1e-12)
{
    check_close(b.xy, xy, what + ".xy", tolerance);
    check_close(b.xz, xz, what + ".xz", tolerance);
    check_close(b.yz, yz, what + ".yz", tolerance);
}

void check_trivector(const Trivector3<double>& t, double e123, const std::string& what,
                     double tolerance = 1e-12)
{
    check_close(t.e123, e123, what + ".e123", tolerance);
}

constexpr double kPi = 3.14159265358979323846;

// Build a multivector from raw components, in grade order.
Multivector3<double> make_mv(double s, double x, double y, double z,
                             double xy, double xz, double yz, double t)
{
    return Multivector3<double>(Scalar<double>(s), Vector3<double>(x, y, z),
                                Bivector3<double>(xy, xz, yz), Trivector3<double>(t));
}

void check_mv(const Multivector3<double>& m, double s, double x, double y, double z,
              double xy, double xz, double yz, double t, const std::string& what)
{
    check_close(m.scalar.value, s, what + ".scalar");
    check_vector(m.vector, x, y, z, what);
    check_bivector(m.bivector, xy, xz, yz, what);
    check_close(m.trivector.e123, t, what + ".trivector");
}

// Total absolute component difference; 0 means identical.
double mv_difference(const Multivector3<double>& a, const Multivector3<double>& b)
{
    return std::fabs(a.scalar.value - b.scalar.value)
         + std::fabs(a.vector.x - b.vector.x)
         + std::fabs(a.vector.y - b.vector.y)
         + std::fabs(a.vector.z - b.vector.z)
         + std::fabs(a.bivector.xy - b.bivector.xy)
         + std::fabs(a.bivector.xz - b.bivector.xz)
         + std::fabs(a.bivector.yz - b.bivector.yz)
         + std::fabs(a.trivector.e123 - b.trivector.e123);
}

Multivector3<double> rotor_to_mv(const Rotor3<double>& r)
{
    return Multivector3<double>(r.scalar, Vector3<double>(), r.bivector, Trivector3<double>());
}

// ---------------------------------------------------------------------------
// Compile-time checks
//
// Everything in the library is constexpr, so these must evaluate during
// compilation. If constexpr-ness regresses, the build fails rather than the run.
// std::sqrt is not portably constexpr, so norm/magnitude stay out of here.
// ---------------------------------------------------------------------------

constexpr Vector3<double> kCtA(1.0, 2.0, 3.0);
constexpr Vector3<double> kCtB(4.0, 5.0, 6.0);
static_assert((kCtA | kCtB).value == 32.0, "dot product must be constexpr");
static_assert((kCtA ^ kCtB).xy == -3.0, "wedge product must be constexpr");
static_assert(Scalar<double>(2.0).value == 2.0, "Scalar construction must be constexpr");
static_assert(Vector3<double>().x == 0.0, "default Vector3 must be constexpr and zeroed");

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

void test_scalar()
{
    section("Scalar");

    check_scalar(Scalar<double>(), 0.0, "default ctor zeroes");
    check_scalar(Scalar<double>(3.5), 3.5, "value ctor");

    const Scalar<double> a(6.0);
    const Scalar<double> b(4.0);

    check_scalar(a + b, 10.0, "operator+");
    check_scalar(a - b, 2.0, "operator-");
    check_scalar(a * b, 24.0, "operator* (Scalar)");
    check_scalar(a / b, 1.5, "operator/");
    check_scalar(-a, -6.0, "unary operator-");

    // Cross-type multiplication: Scalar on the left, a graded type on the right.
    const Scalar<double> two(2.0);
    check_vector(two * Vector3<double>(1.0, 2.0, 3.0), 2.0, 4.0, 6.0, "Scalar * Vector3");
    check_bivector(two * Bivector3<double>(1.0, 2.0, 3.0), 2.0, 4.0, 6.0, "Scalar * Bivector3");
    check_trivector(two * Trivector3<double>(5.0), 10.0, "Scalar * Trivector3");

    // Composite types: the scalar must distribute over every grade.
    const Multivector3<double> mv(
        Scalar<double>(1.0), Vector3<double>(2.0, 3.0, 4.0),
        Bivector3<double>(5.0, 6.0, 7.0), Trivector3<double>(8.0));
    const Multivector3<double> scaledMv = two * mv;
    check_scalar(scaledMv.scalar, 2.0, "Scalar * Multivector3 scalar");
    check_vector(scaledMv.vector, 4.0, 6.0, 8.0, "Scalar * Multivector3 vector");
    check_bivector(scaledMv.bivector, 10.0, 12.0, 14.0, "Scalar * Multivector3 bivector");
    check_trivector(scaledMv.trivector, 16.0, "Scalar * Multivector3 trivector");

    const Rotor3<double> rot(Scalar<double>(1.0), Bivector3<double>(2.0, 3.0, 4.0));
    const Rotor3<double> scaledRot = two * rot;
    check_scalar(scaledRot.scalar, 2.0, "Scalar * Rotor3 scalar");
    check_bivector(scaledRot.bivector, 4.0, 6.0, 8.0, "Scalar * Rotor3 bivector");
}

void test_vector3()
{
    section("Vector3");

    check_vector(Vector3<double>(), 0.0, 0.0, 0.0, "default ctor zeroes");
    check_vector(Vector3<double>(1.0, 2.0, 3.0), 1.0, 2.0, 3.0, "component ctor");

    const Vector3<double> a(1.0, 2.0, 3.0);
    const Vector3<double> b(4.0, 5.0, 6.0);

    check_vector(a + b, 5.0, 7.0, 9.0, "operator+");
    check_vector(a - b, -3.0, -3.0, -3.0, "operator-");
    check_vector(-a, -1.0, -2.0, -3.0, "unary operator-");
    check_vector(a * Scalar<double>(2.0), 2.0, 4.0, 6.0, "operator* (Scalar)");
    check_vector(a / Scalar<double>(2.0), 0.5, 1.0, 1.5, "operator/ (Scalar)");

    // Scalar multiplication agrees across the two spellings.
    const Vector3<double> lhs = Scalar<double>(3.0) * a;
    const Vector3<double> rhs = a * Scalar<double>(3.0);
    check_vector(lhs, rhs.x, rhs.y, rhs.z, "Scalar*Vector3 == Vector3*Scalar");
}

void test_bivector3()
{
    section("Bivector3");

    check_bivector(Bivector3<double>(), 0.0, 0.0, 0.0, "default ctor zeroes");
    check_bivector(Bivector3<double>(1.0, 2.0, 3.0), 1.0, 2.0, 3.0, "component ctor");

    const Bivector3<double> a(1.0, 2.0, 3.0);
    const Bivector3<double> b(4.0, 5.0, 6.0);

    check_bivector(a + b, 5.0, 7.0, 9.0, "operator+");
    check_bivector(a - b, -3.0, -3.0, -3.0, "operator-");
    check_bivector(-a, -1.0, -2.0, -3.0, "unary operator-");
    check_bivector(a * Scalar<double>(2.0), 2.0, 4.0, 6.0, "operator* (Scalar)");
    check_bivector(a / Scalar<double>(2.0), 0.5, 1.0, 1.5, "operator/ (Scalar)");

    check_scalar(a.magnitude(), std::sqrt(14.0), "magnitude");
    check_scalar(Bivector3<double>().magnitude(), 0.0, "magnitude of zero bivector");

    // magnitude() must agree with the free norm() function.
    check_close(a.magnitude().value, CliffordCore::norm(a).value, "magnitude == norm");
}

void test_trivector3()
{
    section("Trivector3");

    check_trivector(Trivector3<double>(), 0.0, "default ctor zeroes");
    check_trivector(Trivector3<double>(7.0), 7.0, "value ctor");

    const Trivector3<double> a(6.0);
    const Trivector3<double> b(3.0);

    check_trivector(a + b, 9.0, "operator+");
    check_trivector(a - b, 3.0, "operator-");
    check_trivector(-a, -6.0, "unary operator-");

    // The geometric product of two pseudoscalars is a SCALAR, not a trivector:
    // e123 squares to -1, so (6 e123)(3 e123) = -18. Division is a * inverse(b),
    // where the two minus signs cancel and leave a plain 6/3.
    static_assert(std::is_same<decltype(a * b), Scalar<double>>::value,
                  "Trivector3 * Trivector3 produces a Scalar");
    static_assert(std::is_same<decltype(a / b), Scalar<double>>::value,
                  "Trivector3 / Trivector3 produces a Scalar");
    check_scalar(a * b, -18.0, "operator* (Trivector3) is -ab");
    check_scalar(a / b, 2.0, "operator/ (Trivector3)");

    // The defining case, and the one the old component-wise version got wrong.
    check_scalar(Trivector3<double>(1.0) * Trivector3<double>(1.0), -1.0, "e123 * e123 == -1");

    // Division must undo multiplication: (a/b) scaled back by b returns a.
    check_trivector((a / b) * b, a.e123, "(a/b) * b == a");

    check_trivector(a * Scalar<double>(2.0), 12.0, "operator* (Scalar)");
    check_trivector(a / Scalar<double>(2.0), 3.0, "operator/ (Scalar)");

    // Raw numeric overloads exist so `t * 2.0` stays unambiguous now that both
    // Trivector3 and Scalar are constructible from one T. They must agree with
    // the Scalar spelling.
    check_trivector(a * 2.0, 12.0, "operator* (raw T)");
    check_trivector(a / 2.0, 3.0, "operator/ (raw T)");
    check_close((a * 2.0).e123, (a * Scalar<double>(2.0)).e123, "raw T == Scalar overload");

    check_scalar(a.magnitude(), 6.0, "magnitude");
    check_scalar(Trivector3<double>(-6.0).magnitude(), 6.0, "magnitude is non-negative");
}

void test_multivector3()
{
    section("Multivector3");

    const Multivector3<double> zero;
    check_scalar(zero.scalar, 0.0, "default ctor scalar");
    check_vector(zero.vector, 0.0, 0.0, 0.0, "default ctor vector");
    check_bivector(zero.bivector, 0.0, 0.0, 0.0, "default ctor bivector");
    check_trivector(zero.trivector, 0.0, "default ctor trivector");

    const Multivector3<double> a(
        Scalar<double>(1.0), Vector3<double>(2.0, 3.0, 4.0),
        Bivector3<double>(5.0, 6.0, 7.0), Trivector3<double>(8.0));
    const Multivector3<double> b(
        Scalar<double>(10.0), Vector3<double>(20.0, 30.0, 40.0),
        Bivector3<double>(50.0, 60.0, 70.0), Trivector3<double>(80.0));

    const Multivector3<double> sum = a + b;
    check_scalar(sum.scalar, 11.0, "operator+ scalar");
    check_vector(sum.vector, 22.0, 33.0, 44.0, "operator+ vector");
    check_bivector(sum.bivector, 55.0, 66.0, 77.0, "operator+ bivector");
    check_trivector(sum.trivector, 88.0, "operator+ trivector");

    const Multivector3<double> diff = b - a;
    check_scalar(diff.scalar, 9.0, "operator- scalar");
    check_vector(diff.vector, 18.0, 27.0, 36.0, "operator- vector");
    check_bivector(diff.bivector, 45.0, 54.0, 63.0, "operator- bivector");
    check_trivector(diff.trivector, 72.0, "operator- trivector");

    const Multivector3<double> negated = -a;
    check_scalar(negated.scalar, -1.0, "unary operator- scalar");
    check_vector(negated.vector, -2.0, -3.0, -4.0, "unary operator- vector");
    check_bivector(negated.bivector, -5.0, -6.0, -7.0, "unary operator- bivector");
    check_trivector(negated.trivector, -8.0, "unary operator- trivector");

    const Multivector3<double> scaled = a * Scalar<double>(2.0);
    check_scalar(scaled.scalar, 2.0, "operator* (Scalar) scalar");
    check_vector(scaled.vector, 4.0, 6.0, 8.0, "operator* (Scalar) vector");
    check_bivector(scaled.bivector, 10.0, 12.0, 14.0, "operator* (Scalar) bivector");
    check_trivector(scaled.trivector, 16.0, "operator* (Scalar) trivector");

    const Multivector3<double> halved = a / Scalar<double>(2.0);
    check_scalar(halved.scalar, 0.5, "operator/ (Scalar) scalar");
    check_vector(halved.vector, 1.0, 1.5, 2.0, "operator/ (Scalar) vector");
    check_bivector(halved.bivector, 2.5, 3.0, 3.5, "operator/ (Scalar) bivector");
    check_trivector(halved.trivector, 4.0, "operator/ (Scalar) trivector");

    // Scaling agrees whichever side the scalar sits on.
    const Multivector3<double> fromLeft = Scalar<double>(2.0) * a;
    check_close(fromLeft.scalar.value, scaled.scalar.value, "Scalar*MV == MV*Scalar scalar");
    check_vector(fromLeft.vector, scaled.vector.x, scaled.vector.y, scaled.vector.z,
                 "Scalar*MV == MV*Scalar vector");
    check_trivector(fromLeft.trivector, scaled.trivector.e123, "Scalar*MV == MV*Scalar trivector");

    // Negation is subtraction from zero.
    const Multivector3<double> viaSubtraction = Multivector3<double>() - a;
    check_scalar(viaSubtraction.scalar, negated.scalar.value, "-a == 0 - a scalar");
    check_trivector(viaSubtraction.trivector, negated.trivector.e123, "-a == 0 - a trivector");
}

void test_rotor3()
{
    section("Rotor3");

    const Rotor3<double> zero;
    check_scalar(zero.scalar, 0.0, "default ctor scalar");
    check_bivector(zero.bivector, 0.0, 0.0, 0.0, "default ctor bivector");

    const Rotor3<double> a(Scalar<double>(1.0), Bivector3<double>(2.0, 3.0, 4.0));
    const Rotor3<double> b(Scalar<double>(10.0), Bivector3<double>(20.0, 30.0, 40.0));

    const Rotor3<double> sum = a + b;
    check_scalar(sum.scalar, 11.0, "operator+ scalar");
    check_bivector(sum.bivector, 22.0, 33.0, 44.0, "operator+ bivector");

    const Rotor3<double> diff = b - a;
    check_scalar(diff.scalar, 9.0, "operator- scalar");
    check_bivector(diff.bivector, 18.0, 27.0, 36.0, "operator- bivector");

    const Rotor3<double> scaled = a * Scalar<double>(2.0);
    check_scalar(scaled.scalar, 2.0, "operator* (Scalar) scalar");
    check_bivector(scaled.bivector, 4.0, 6.0, 8.0, "operator* (Scalar) bivector");
}

// ---------------------------------------------------------------------------
// Operations
// ---------------------------------------------------------------------------

void test_dot_product()
{
    section("dot_product");

    const Vector3<double> a(1.0, 2.0, 3.0);
    const Vector3<double> b(4.0, 5.0, 6.0);

    check_scalar(CliffordCore::dot_product(a, b), 32.0, "dot_product");
    check_scalar(a | b, 32.0, "operator|");
    check_close((a | b).value, CliffordCore::dot_product(a, b).value, "operator| == dot_product");

    // Algebraic properties, not just one hardcoded number.
    check_close((a | b).value, (b | a).value, "dot product is commutative");
    check_scalar(a | a, 14.0, "v|v equals squared norm");

    const Vector3<double> x(1.0, 0.0, 0.0);
    const Vector3<double> y(0.0, 1.0, 0.0);
    check_scalar(x | y, 0.0, "orthogonal basis vectors have zero dot");
}

void test_wedge_product()
{
    section("wedge_product");

    const Vector3<double> a(1.0, 2.0, 3.0);
    const Vector3<double> b(4.0, 5.0, 6.0);

    check_bivector(CliffordCore::wedge_product(a, b), -3.0, -6.0, -3.0, "wedge_product");
    check_bivector(a ^ b, -3.0, -6.0, -3.0, "operator^");

    // Antisymmetry: a^b == -(b^a).
    const Bivector3<double> ab = a ^ b;
    const Bivector3<double> ba = b ^ a;
    check_bivector(ba, -ab.xy, -ab.xz, -ab.yz, "wedge is antisymmetric");

    // A vector wedged with itself is zero. Each component is x*y - y*x with
    // identical operands, so this is exactly zero and can be compared with ==.
    const Bivector3<double> selfWedge = a ^ a;
    check(selfWedge.xy == 0.0 && selfWedge.xz == 0.0 && selfWedge.yz == 0.0,
          "v^v is exactly zero");

    // Basis check, pinning the (xy, xz, yz) component convention.
    const Vector3<double> x(1.0, 0.0, 0.0);
    const Vector3<double> y(0.0, 1.0, 0.0);
    check_bivector(x ^ y, 1.0, 0.0, 0.0, "ex^ey is the xy basis bivector");
}

void test_geometric_product()
{
    section("geometric_product");

    const Vector3<double> a(1.0, 2.0, 3.0);
    const Vector3<double> b(4.0, 5.0, 6.0);

    const Multivector3<double> gp = CliffordCore::geometric_product(a, b);
    check_scalar(gp.scalar, 32.0, "scalar part is the dot product");
    check_bivector(gp.bivector, -3.0, -6.0, -3.0, "bivector part is the wedge product");
    check_vector(gp.vector, 0.0, 0.0, 0.0, "vector part is zero");
    check_trivector(gp.trivector, 0.0, "trivector part is zero");

    // operator* must agree with the named function.
    const Multivector3<double> viaOperator = a * b;
    check_close(viaOperator.scalar.value, gp.scalar.value, "operator* scalar matches");
    check_bivector(viaOperator.bivector, gp.bivector.xy, gp.bivector.xz, gp.bivector.yz,
                   "operator* bivector matches");

    // A vector times itself is a pure scalar equal to its squared norm.
    const Multivector3<double> square = a * a;
    check_scalar(square.scalar, 14.0, "v*v scalar part is squared norm");
    check_bivector(square.bivector, 0.0, 0.0, 0.0, "v*v has no bivector part");

    // rotor_product packs the same numbers into a Rotor3.
    const Rotor3<double> rotor = CliffordCore::rotor_product(a, b);
    check_scalar(rotor.scalar, 32.0, "rotor_product scalar");
    check_bivector(rotor.bivector, -3.0, -6.0, -3.0, "rotor_product bivector");
}

void test_norm()
{
    section("norm");

    const Vector3<double> a(1.0, 2.0, 3.0);
    check_scalar(CliffordCore::norm(a), std::sqrt(14.0), "norm(Vector3)");
    check_scalar(CliffordCore::squared_norm(a), 14.0, "squared_norm(Vector3)");

    const Bivector3<double> b(1.0, 2.0, 3.0);
    check_scalar(CliffordCore::norm(b), std::sqrt(14.0), "norm(Bivector3)");
    check_scalar(CliffordCore::squared_norm(b), 14.0, "squared_norm(Bivector3)");

    // norm squared must equal squared_norm.
    check_close(CliffordCore::norm(a).value * CliffordCore::norm(a).value,
                CliffordCore::squared_norm(a).value, "norm^2 == squared_norm");

    check_scalar(CliffordCore::norm(Vector3<double>()), 0.0, "norm of zero vector");
    check_scalar(CliffordCore::norm(Vector3<double>(3.0, 4.0, 0.0)), 5.0, "3-4-5 triangle");
}

void test_inverse()
{
    section("inverse");

    const Vector3<double> a(1.0, 2.0, 3.0);
    const Vector3<double> inv = CliffordCore::inverse(a);
    check_vector(inv, 1.0 / 14.0, 2.0 / 14.0, 3.0 / 14.0, "inverse(Vector3)");

    // Defining property: v * v^-1 == 1, a pure scalar of value one.
    const Multivector3<double> product = a * inv;
    check_scalar(product.scalar, 1.0, "v * inverse(v) has unit scalar part");
    check_bivector(product.bivector, 0.0, 0.0, 0.0, "v * inverse(v) has no bivector part");

    // A unit vector is its own inverse.
    const Vector3<double> unit(1.0, 0.0, 0.0);
    check_vector(CliffordCore::inverse(unit), 1.0, 0.0, 0.0, "unit vector is its own inverse");
}

void test_reverse()
{
    section("reverse");

    const Vector3<double> v(1.0, 2.0, 3.0);
    check_vector(CliffordCore::reverse(v), 1.0, 2.0, 3.0, "reverse(Vector3) is identity");

    const Bivector3<double> b(1.0, 2.0, 3.0);
    check_bivector(CliffordCore::reverse(b), -1.0, -2.0, -3.0, "reverse(Bivector3) negates");

    // Reverse is an involution: applying it twice returns the original.
    const Bivector3<double> twice = CliffordCore::reverse(CliffordCore::reverse(b));
    check_bivector(twice, b.xy, b.xz, b.yz, "reverse(reverse(b)) == b");
}

void test_dual()
{
    section("dual");

    // dual(a) is a * e123. Checked here against the literal basis mapping, and
    // below against the general multivector product.
    const Vector3<double> v(1.0, 2.0, 3.0);
    check_bivector(CliffordCore::dual(v), 3.0, -2.0, 1.0, "dual(Vector3) maps x,y,z to yz,-xz,xy");

    const Bivector3<double> b(4.0, 5.0, 6.0);
    check_vector(CliffordCore::dual(b), -6.0, 5.0, -4.0, "dual(Bivector3)");

    // Grade 0 and grade 3 are dual to each other too.
    check_trivector(CliffordCore::dual(Scalar<double>(7.0)), 7.0, "dual(Scalar) is the pseudoscalar");
    check_scalar(CliffordCore::dual(Trivector3<double>(7.0)), -7.0, "dual(Trivector3) negates");

    // Basis checks, pinning each individual mapping.
    check_bivector(CliffordCore::dual(Vector3<double>(1, 0, 0)), 0.0, 0.0, 1.0, "dual(e1) = e23");
    check_bivector(CliffordCore::dual(Vector3<double>(0, 1, 0)), 0.0, -1.0, 0.0, "dual(e2) = -e13");
    check_bivector(CliffordCore::dual(Vector3<double>(0, 0, 1)), 1.0, 0.0, 0.0, "dual(e3) = e12");

    // The pseudoscalar squares to -1, so the dual is an involution up to sign:
    // applying it twice negates rather than returning the original.
    const Vector3<double> twice = CliffordCore::dual(CliffordCore::dual(v));
    check_vector(twice, -v.x, -v.y, -v.z, "dual(dual(v)) == -v");
    check_scalar(CliffordCore::dual(CliffordCore::dual(Scalar<double>(7.0))), -7.0,
                 "dual(dual(s)) == -s");

    // The definition itself: dual(a) must equal a * e123 under the general
    // product. This is what makes it the algebraic dual rather than a relabelling.
    const Multivector3<double> I = make_mv(0, 0, 0, 0, 0, 0, 0, 1);
    const Multivector3<double> vAsMv = make_mv(0, v.x, v.y, v.z, 0, 0, 0, 0);
    const Bivector3<double> dualV = CliffordCore::dual(v);
    check_close(mv_difference(vAsMv * I, make_mv(0, 0, 0, 0, dualV.xy, dualV.xz, dualV.yz, 0)), 0.0,
                "dual(v) == v * e123");

    const Multivector3<double> bAsMv = make_mv(0, 0, 0, 0, b.xy, b.xz, b.yz, 0);
    const Vector3<double> dualB = CliffordCore::dual(b);
    check_close(mv_difference(bAsMv * I, make_mv(0, dualB.x, dualB.y, dualB.z, 0, 0, 0, 0)), 0.0,
                "dual(b) == b * e123");

    // The dual of a wedge is the cross product, which is the classic reason to
    // want a dual in 3D: e1 ^ e2 is the xy plane, whose dual is the e3 axis.
    const Vector3<double> e1(1, 0, 0);
    const Vector3<double> e2(0, 1, 0);
    check_vector(CliffordCore::dual(e1 ^ e2), 0.0, 0.0, -1.0, "dual(e1 ^ e2) is along e3");
}

void test_exp_log()
{
    section("exp/log");

    // exp of the zero bivector is the identity rotor.
    const Rotor3<double> identity = CliffordCore::exp(Bivector3<double>(0.0, 0.0, 0.0));
    check_scalar(identity.scalar, 1.0, "exp(0) scalar");
    check_bivector(identity.bivector, 0.0, 0.0, 0.0, "exp(0) bivector");

    // A quarter turn in the xy plane: cos(pi/2) + sin(pi/2) * xy.
    const Rotor3<double> quarter = CliffordCore::exp(Bivector3<double>(kPi / 2.0, 0.0, 0.0));
    check_close(quarter.scalar.value, 0.0, "exp(pi/2 xy) scalar", 1e-15);
    check_bivector(quarter.bivector, 1.0, 0.0, 0.0, "exp(pi/2 xy) bivector");

    // exp must always produce a unit rotor: scalar^2 + |bivector|^2 == 1.
    const Bivector3<double> b(0.3, 0.4, 0.5);
    const Rotor3<double> r = CliffordCore::exp(b);
    check_close(r.scalar.value * r.scalar.value + CliffordCore::squared_norm(r.bivector).value,
                1.0, "exp produces a unit rotor", 1e-15);

    // log of the identity rotor is the zero bivector.
    const Bivector3<double> logIdentity =
        CliffordCore::log(Rotor3<double>(Scalar<double>(1.0), Bivector3<double>(0.0, 0.0, 0.0)));
    check_bivector(logIdentity, 0.0, 0.0, 0.0, "log(identity) is zero");

    // log and exp are inverses. Trig costs precision, so allow a looser bound.
    const Bivector3<double> there = CliffordCore::log(CliffordCore::exp(b));
    check_bivector(there, b.xy, b.xz, b.yz, "log(exp(b)) == b");

    const Rotor3<double> back = CliffordCore::exp(CliffordCore::log(r));
    check_close(back.scalar.value, r.scalar.value, "exp(log(r)) scalar", 1e-14);
    check_bivector(back.bivector, r.bivector.xy, r.bivector.xz, r.bivector.yz,
                   "exp(log(r)) bivector");

    // A rotor whose scalar has drifted just outside [-1, 1] must not yield NaN:
    // std::acos would be out of domain without the clamp in log().
    const Bivector3<double> drifted =
        CliffordCore::log(Rotor3<double>(Scalar<double>(1.0 + 2e-16), Bivector3<double>(1.0, 0.0, 0.0)));
    check(!std::isnan(drifted.xy) && !std::isnan(drifted.xz) && !std::isnan(drifted.yz),
          "log clamps scalar drift above 1 instead of returning NaN");

    const Bivector3<double> driftedLow =
        CliffordCore::log(Rotor3<double>(Scalar<double>(-1.0 - 2e-16), Bivector3<double>(1.0, 0.0, 0.0)));
    check(!std::isnan(driftedLow.xy), "log clamps scalar drift below -1");
}

// ---------------------------------------------------------------------------
// The algebra itself
// ---------------------------------------------------------------------------

void test_basis_multiplication_table()
{
    section("basis table");

    // The whole algebra follows from these. Basis: e_i^2 = +1, bivectors stored
    // as (xy, xz, yz) = (e1e2, e1e3, e2e3), pseudoscalar e123 = e1e2e3.
    const Multivector3<double> e1  = make_mv(0, 1, 0, 0, 0, 0, 0, 0);
    const Multivector3<double> e2  = make_mv(0, 0, 1, 0, 0, 0, 0, 0);
    const Multivector3<double> e3  = make_mv(0, 0, 0, 1, 0, 0, 0, 0);
    const Multivector3<double> e12 = make_mv(0, 0, 0, 0, 1, 0, 0, 0);
    const Multivector3<double> e13 = make_mv(0, 0, 0, 0, 0, 1, 0, 0);
    const Multivector3<double> e23 = make_mv(0, 0, 0, 0, 0, 0, 1, 0);
    const Multivector3<double> I   = make_mv(0, 0, 0, 0, 0, 0, 0, 1);

    // Vectors square to +1.
    check_close((e1 * e1).scalar.value, 1.0, "e1^2 = +1");
    check_close((e2 * e2).scalar.value, 1.0, "e2^2 = +1");
    check_close((e3 * e3).scalar.value, 1.0, "e3^2 = +1");

    // Bivectors and the pseudoscalar square to -1.
    check_close((e12 * e12).scalar.value, -1.0, "e12^2 = -1");
    check_close((e13 * e13).scalar.value, -1.0, "e13^2 = -1");
    check_close((e23 * e23).scalar.value, -1.0, "e23^2 = -1");
    check_close((I * I).scalar.value, -1.0, "e123^2 = -1");

    // Distinct basis vectors anticommute.
    check_mv(e1 * e2, 0, 0, 0, 0, 1, 0, 0, 0, "e1e2 = e12");
    check_mv(e2 * e1, 0, 0, 0, 0, -1, 0, 0, 0, "e2e1 = -e12");
    check_mv(e1 * e3, 0, 0, 0, 0, 0, 1, 0, 0, "e1e3 = e13");
    check_mv(e2 * e3, 0, 0, 0, 0, 0, 0, 1, 0, "e2e3 = e23");
    check_mv(e1 * e2 * e3, 0, 0, 0, 0, 0, 0, 0, 1, "e1e2e3 = e123");

    // Bivector cross products. These signs are what the (xy, xz, yz) ordering
    // buys, and they are the easiest thing in the library to get backwards.
    check_mv(e12 * e13, 0, 0, 0, 0, 0, 0, -1, 0, "e12 e13 = -e23");
    check_mv(e13 * e12, 0, 0, 0, 0, 0, 0, 1, 0, "e13 e12 = +e23");
    check_mv(e12 * e23, 0, 0, 0, 0, 0, 1, 0, 0, "e12 e23 = +e13");
    check_mv(e23 * e12, 0, 0, 0, 0, 0, -1, 0, 0, "e23 e12 = -e13");
    check_mv(e13 * e23, 0, 0, 0, 0, -1, 0, 0, 0, "e13 e23 = -e12");
    check_mv(e23 * e13, 0, 0, 0, 0, 1, 0, 0, 0, "e23 e13 = +e12");

    // The pseudoscalar is central in 3D: it commutes with everything.
    check_close(mv_difference(e1 * I, I * e1), 0.0, "e123 commutes with e1");
    check_close(mv_difference(e12 * I, I * e12), 0.0, "e123 commutes with e12");
}

void test_multivector_product()
{
    section("multivector product");

    const Multivector3<double> a = make_mv(1, 2, -3, 4, -5, 6, 7, -8);
    const Multivector3<double> b = make_mv(-2, 1, 5, -3, 2, -4, 1, 6);
    const Multivector3<double> c = make_mv(3, -1, 2, 1, -2, 3, -1, 2);

    // Associativity is the strongest single check on the multiplication table:
    // a wrong sign anywhere almost certainly breaks it.
    check_close(mv_difference((a * b) * c, a * (b * c)), 0.0,
                "product is associative", 1e-10);

    // It is NOT commutative -- that is the point of a geometric algebra.
    check(mv_difference(a * b, b * a) > 1e-9, "product is not commutative");

    // Distributes over addition.
    check_close(mv_difference(a * (b + c), (a * b) + (a * c)), 0.0,
                "product distributes over addition", 1e-10);

    // The scalar 1 is the identity.
    const Multivector3<double> one = make_mv(1, 0, 0, 0, 0, 0, 0, 0);
    check_close(mv_difference(one * a, a), 0.0, "1 * a == a");
    check_close(mv_difference(a * one, a), 0.0, "a * 1 == a");

    // The general product must reproduce the specialised vector-vector one.
    const Vector3<double> u(1, 2, 3);
    const Vector3<double> w(4, 5, 6);
    const Multivector3<double> viaGeneral =
        make_mv(0, 1, 2, 3, 0, 0, 0, 0) * make_mv(0, 4, 5, 6, 0, 0, 0, 0);
    check_close(mv_difference(viaGeneral, CliffordCore::geometric_product(u, w)), 0.0,
                "general product matches geometric_product(v, v)");

    // ...and the specialised rotor-rotor one.
    const Rotor3<double> r1(Scalar<double>(0.5), Bivector3<double>(0.3, -0.2, 0.7));
    const Rotor3<double> r2(Scalar<double>(-0.4), Bivector3<double>(0.1, 0.6, -0.5));
    check_close(mv_difference(rotor_to_mv(r1) * rotor_to_mv(r2), rotor_to_mv(r1 * r2)), 0.0,
                "general product matches Rotor3 * Rotor3");

    // A vector times itself is its squared norm, with nothing else left over.
    check_close(mv_difference(make_mv(0, 1, 2, 3, 0, 0, 0, 0) * make_mv(0, 1, 2, 3, 0, 0, 0, 0),
                              make_mv(14, 0, 0, 0, 0, 0, 0, 0)), 0.0, "v * v == |v|^2");
}

void test_norm_all_grades()
{
    section("norm (all grades)");

    check_scalar(CliffordCore::norm(Scalar<double>(-3.0)), 3.0, "norm(Scalar) is absolute");
    check_scalar(CliffordCore::norm(Trivector3<double>(-4.0)), 4.0, "norm(Trivector3)");
    check_scalar(CliffordCore::squared_norm(Trivector3<double>(4.0)), 16.0, "squared_norm(Trivector3)");

    // 1 + 4 + 9 + 16 + 25 + 36 + 49 + 64 = 204
    const Multivector3<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8);
    check_scalar(CliffordCore::squared_norm(m), 204.0, "squared_norm(Multivector3)");
    check_scalar(CliffordCore::norm(m), std::sqrt(204.0), "norm(Multivector3)");

    const Rotor3<double> r(Scalar<double>(3.0), Bivector3<double>(4.0, 0.0, 0.0));
    check_scalar(CliffordCore::squared_norm(r), 25.0, "squared_norm(Rotor3)");
    check_scalar(CliffordCore::norm(r), 5.0, "norm(Rotor3)");

    check_scalar(CliffordCore::norm(Multivector3<double>()), 0.0, "norm of zero multivector");
}

void test_reverse_all_grades()
{
    section("reverse (all grades)");

    // Reverse flips grade k by (-1)^(k(k-1)/2): grades 0 and 1 keep their sign,
    // grades 2 and 3 are negated.
    check_trivector(CliffordCore::reverse(Trivector3<double>(5.0)), -5.0, "reverse(Trivector3) negates");

    const Multivector3<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8);
    check_mv(CliffordCore::reverse(m), 1, 2, 3, 4, -5, -6, -7, -8, "reverse(Multivector3)");

    // Reverse is an involution on every type.
    check_close(mv_difference(CliffordCore::reverse(CliffordCore::reverse(m)), m), 0.0,
                "reverse(reverse(m)) == m");
    check_trivector(CliffordCore::reverse(CliffordCore::reverse(Trivector3<double>(5.0))), 5.0,
                    "reverse(reverse(t)) == t");

    // reverse(ab) == reverse(b) reverse(a) -- the defining anti-automorphism.
    const Multivector3<double> a = make_mv(1, 2, -3, 4, -5, 6, 7, -8);
    const Multivector3<double> b = make_mv(-2, 1, 5, -3, 2, -4, 1, 6);
    check_close(mv_difference(CliffordCore::reverse(a * b),
                              CliffordCore::reverse(b) * CliffordCore::reverse(a)), 0.0,
                "reverse(ab) == reverse(b) reverse(a)", 1e-10);
}

void test_inverse_all_grades()
{
    section("inverse (all grades)");

    check_scalar(CliffordCore::inverse(Scalar<double>(4.0)), 0.25, "inverse(Scalar)");

    // Each inverse is checked by its defining property, x * inverse(x) == 1,
    // evaluated with the general product rather than by restating the formula.
    const Multivector3<double> one = make_mv(1, 0, 0, 0, 0, 0, 0, 0);

    const Vector3<double> v(1.0, 2.0, 3.0);
    check_close(mv_difference(make_mv(0, v.x, v.y, v.z, 0, 0, 0, 0)
                    * make_mv(0, CliffordCore::inverse(v).x, CliffordCore::inverse(v).y,
                              CliffordCore::inverse(v).z, 0, 0, 0, 0), one),
                0.0, "v * inverse(v) == 1");

    // A bivector squares to -|b|^2, so its inverse carries a minus sign.
    const Bivector3<double> b(1.0, 2.0, 3.0);
    const Bivector3<double> bi = CliffordCore::inverse(b);
    check_bivector(bi, -1.0 / 14.0, -2.0 / 14.0, -3.0 / 14.0, "inverse(Bivector3) is -b/|b|^2");
    check_close(mv_difference(make_mv(0, 0, 0, 0, b.xy, b.xz, b.yz, 0)
                    * make_mv(0, 0, 0, 0, bi.xy, bi.xz, bi.yz, 0), one),
                0.0, "b * inverse(b) == 1");

    // The pseudoscalar also squares to -1, so the same minus sign applies.
    const Trivector3<double> t(4.0);
    const Trivector3<double> ti = CliffordCore::inverse(t);
    check_trivector(ti, -0.25, "inverse(Trivector3) is -t/|t|^2");
    check_close(mv_difference(make_mv(0, 0, 0, 0, 0, 0, 0, t.e123)
                    * make_mv(0, 0, 0, 0, 0, 0, 0, ti.e123), one),
                0.0, "t * inverse(t) == 1");

    // A rotor inverts by reversing, not by dividing itself by its norm.
    const Rotor3<double> r(Scalar<double>(0.6), Bivector3<double>(0.8, 0.0, 0.0));
    const Rotor3<double> ri = CliffordCore::inverse(r);
    const Rotor3<double> product = r * ri;
    check_scalar(product.scalar, 1.0, "r * inverse(r) scalar is 1");
    check_bivector(product.bivector, 0.0, 0.0, 0.0, "r * inverse(r) has no bivector part");
    // For a unit rotor the inverse is exactly the reverse.
    check_bivector(ri.bivector, CliffordCore::reverse(r).bivector.xy,
                   CliffordCore::reverse(r).bivector.xz,
                   CliffordCore::reverse(r).bivector.yz,
                   "unit rotor inverse == reverse");

    // The general multivector inverse. m/|m|^2 is not it; the conjugate route is.
    const Multivector3<double> m = make_mv(1, 2, -3, 4, -5, 6, 7, -8);
    check_close(mv_difference(m * CliffordCore::inverse(m), one), 0.0,
                "m * inverse(m) == 1", 1e-10);
    check_close(mv_difference(CliffordCore::inverse(m) * m, one), 0.0,
                "inverse(m) * m == 1", 1e-10);

    // Clifford conjugation negates grades 1 and 2 only.
    check_mv(CliffordCore::conjugate(make_mv(1, 2, 3, 4, 5, 6, 7, 8)),
             1, -2, -3, -4, -5, -6, -7, 8, "conjugate(m)");
}

void test_normalize()
{
    section("normalize");

    check_vector(CliffordCore::normalize(Vector3<double>(3.0, 4.0, 0.0)), 0.6, 0.8, 0.0,
                 "normalize(Vector3)");
    check_scalar(CliffordCore::norm(CliffordCore::normalize(Vector3<double>(1.0, 2.0, 3.0))), 1.0,
                 "normalized vector has unit norm");
    check_scalar(CliffordCore::norm(CliffordCore::normalize(Bivector3<double>(1.0, 2.0, 3.0))), 1.0,
                 "normalized bivector has unit norm");
    check_trivector(CliffordCore::normalize(Trivector3<double>(-4.0)), -1.0, "normalize(Trivector3)");
    check_scalar(CliffordCore::norm(CliffordCore::normalize(make_mv(1, 2, 3, 4, 5, 6, 7, 8))), 1.0,
                 "normalized multivector has unit norm");

    // Rotors must be unit length to represent a rotation.
    const Rotor3<double> drifted(Scalar<double>(3.0), Bivector3<double>(4.0, 0.0, 0.0));
    check_scalar(CliffordCore::norm(CliffordCore::normalize(drifted)), 1.0,
                 "normalized rotor has unit norm");
    check_scalar(CliffordCore::normalize(drifted).scalar, 0.6, "normalize(Rotor3) scalar");

    // Zero has no direction to preserve, so it is returned unchanged rather than
    // producing NaN.
    check_vector(CliffordCore::normalize(Vector3<double>()), 0.0, 0.0, 0.0,
                 "normalize(zero vector) stays zero");
    check_scalar(CliffordCore::normalize(Rotor3<double>()).scalar, 0.0,
                 "normalize(zero rotor) stays zero");
}

void test_rotation()
{
    section("rotation");

    // exp(theta * B) rotates by 2*theta in the plane B, so pi/4 gives a quarter
    // turn. e12 spans the xy plane, so it turns about z and must leave e3 alone.
    const Rotor3<double> quarterXY = CliffordCore::exp(Bivector3<double>(kPi / 4.0, 0.0, 0.0));
    check_vector(CliffordCore::rotate(Vector3<double>(1, 0, 0), quarterXY), 0.0, -1.0, 0.0,
                 "e1 rotates to -e2 in the xy plane", 1e-15);
    check_vector(CliffordCore::rotate(Vector3<double>(0, 0, 1), quarterXY), 0.0, 0.0, 1.0,
                 "e3 is the axis of an xy rotation and is fixed", 1e-15);

    // e23 spans the yz plane, so it turns about x and must leave e1 alone.
    const Rotor3<double> quarterYZ = CliffordCore::exp(Bivector3<double>(0.0, 0.0, kPi / 4.0));
    check_vector(CliffordCore::rotate(Vector3<double>(1, 0, 0), quarterYZ), 1.0, 0.0, 0.0,
                 "e1 is the axis of a yz rotation and is fixed", 1e-15);

    // Rotation preserves length.
    const Vector3<double> v(1.0, 2.0, 3.0);
    check_scalar(CliffordCore::norm(CliffordCore::rotate(v, quarterXY)), std::sqrt(14.0),
                 "rotation preserves length");

    // The identity rotor leaves everything alone.
    const Rotor3<double> identity(Scalar<double>(1.0), Bivector3<double>(0.0, 0.0, 0.0));
    check_vector(CliffordCore::rotate(v, identity), v.x, v.y, v.z, "identity rotor fixes v");

    // The closed form must agree with the general product R v reverse(R).
    const Rotor3<double> r =
        CliffordCore::normalize(Rotor3<double>(Scalar<double>(0.3), Bivector3<double>(0.5, -0.7, 0.2)));
    const Multivector3<double> viaGeneral =
        rotor_to_mv(r) * make_mv(0, v.x, v.y, v.z, 0, 0, 0, 0) * rotor_to_mv(CliffordCore::reverse(r));
    const Vector3<double> viaClosedForm = CliffordCore::sandwich(v, r);
    check_vector(viaClosedForm, viaGeneral.vector.x, viaGeneral.vector.y, viaGeneral.vector.z,
                 "sandwich matches R v reverse(R)", 1e-14);
    // The general product must come out a pure vector.
    check_close(viaGeneral.scalar.value, 0.0, "sandwich leaves no scalar part", 1e-14);
    check_close(viaGeneral.trivector.e123, 0.0, "sandwich leaves no trivector part", 1e-14);

    // Rotating twice by R equals rotating once by R*R.
    const Vector3<double> twice = CliffordCore::rotate(CliffordCore::rotate(v, quarterXY), quarterXY);
    const Vector3<double> once = CliffordCore::rotate(v, quarterXY * quarterXY);
    check_vector(twice, once.x, once.y, once.z, "rotating twice == rotating by the composed rotor", 1e-14);

    // A rotor and its inverse undo each other.
    const Vector3<double> roundTrip =
        CliffordCore::rotate(CliffordCore::rotate(v, quarterXY), CliffordCore::inverse(quarterXY));
    check_vector(roundTrip, v.x, v.y, v.z, "rotate then rotate by the inverse returns v", 1e-14);
}

void test_mixed_grade_addition()
{
    section("mixed-grade addition");

    const Scalar<double> s(1.0);
    const Vector3<double> v1(1, 0, 0), v2(0, 1, 0), v3(0, 0, 1), v4(1, 1, 1);
    const Bivector3<double> b1(1, 0, 0), b2(0, 1, 0), b3(0, 0, 1);
    const Trivector3<double> t(5.0);

    // A chain across every grade, in one expression.
    const Multivector3<double> m = s + v1 + v2 + v3 + v4 + b1 + b2 + b3 + t;
    check_mv(m, 1, 2, 2, 2, 1, 1, 1, 5, "nine-term mixed chain");

    // Addition is commutative, so the order of the chain must not matter.
    const Multivector3<double> reordered = t + b3 + b2 + b1 + v4 + v3 + v2 + v1 + s;
    check_close(mv_difference(m, reordered), 0.0, "chain order does not matter");

    // Every mixed pair lands in the right slot.
    check_mv(s + v1, 1, 1, 0, 0, 0, 0, 0, 0, "Scalar + Vector3");
    check_mv(v1 + s, 1, 1, 0, 0, 0, 0, 0, 0, "Vector3 + Scalar");
    check_mv(s + b1, 1, 0, 0, 0, 1, 0, 0, 0, "Scalar + Bivector3");
    check_mv(s + t, 1, 0, 0, 0, 0, 0, 0, 5, "Scalar + Trivector3");
    check_mv(v1 + b1, 0, 1, 0, 0, 1, 0, 0, 0, "Vector3 + Bivector3");
    check_mv(v1 + t, 0, 1, 0, 0, 0, 0, 0, 5, "Vector3 + Trivector3");
    check_mv(b1 + t, 0, 0, 0, 0, 1, 0, 0, 5, "Bivector3 + Trivector3");
    check_mv(m + v1, 1, 3, 2, 2, 1, 1, 1, 5, "Multivector3 + Vector3");
    check_mv(v1 + m, 1, 3, 2, 2, 1, 1, 1, 5, "Vector3 + Multivector3");

    // Same-grade sums must not widen.
    static_assert(std::is_same<decltype(v1 + v2), Vector3<double>>::value,
                  "Vector3 + Vector3 stays a Vector3");
    static_assert(std::is_same<decltype(s + b1), Multivector3<double>>::value,
                  "Scalar + Bivector3 widens to Multivector3");

    // rotor_sum packs the same numbers as s + b, and to_rotor narrows.
    const Rotor3<double> viaSum = CliffordCore::rotor_sum(s, b1);
    check_scalar(viaSum.scalar, 1.0, "rotor_sum scalar");
    check_bivector(viaSum.bivector, 1.0, 0.0, 0.0, "rotor_sum bivector");
    const Rotor3<double> viaNarrow = CliffordCore::to_rotor(s + b1);
    check_scalar(viaNarrow.scalar, 1.0, "to_rotor scalar");
    check_bivector(viaNarrow.bivector, 1.0, 0.0, 0.0, "to_rotor bivector");
    check_mv(CliffordCore::to_multivector(viaSum), 1, 0, 0, 0, 1, 0, 0, 0, "to_multivector(rotor)");
}

void test_mixed_grade_subtraction()
{
    section("mixed-grade subtraction");

    const Scalar<double> s(2.0);
    const Vector3<double> v(1, 2, 3);
    const Bivector3<double> b(4, 5, 6);
    const Trivector3<double> t(7.0);
    const Multivector3<double> m = make_mv(1, 1, 2, 3, 1, 2, 3, 4);

    check_mv(s - v, 2, -1, -2, -3, 0, 0, 0, 0, "Scalar - Vector3");
    check_mv(v - s, -2, 1, 2, 3, 0, 0, 0, 0, "Vector3 - Scalar");
    check_mv(s - b, 2, 0, 0, 0, -4, -5, -6, 0, "Scalar - Bivector3");
    check_mv(s - t, 2, 0, 0, 0, 0, 0, 0, -7, "Scalar - Trivector3");
    check_mv(v - b, 0, 1, 2, 3, -4, -5, -6, 0, "Vector3 - Bivector3");
    check_mv(b - v, 0, -1, -2, -3, 4, 5, 6, 0, "Bivector3 - Vector3");
    check_mv(v - t, 0, 1, 2, 3, 0, 0, 0, -7, "Vector3 - Trivector3");
    check_mv(b - t, 0, 0, 0, 0, 4, 5, 6, -7, "Bivector3 - Trivector3");
    check_mv(m - v, 1, 0, 0, 0, 1, 2, 3, 4, "Multivector3 - Vector3");
    check_mv(v - m, -1, 0, 0, 0, -1, -2, -3, -4, "Vector3 - Multivector3");

    // a - b must equal a + (-b), and a - a must vanish.
    check_close(mv_difference(v - b, v + (-b)), 0.0, "a - b == a + (-b)");
    check_close(mv_difference(m - m, Multivector3<double>()), 0.0, "m - m == 0");

    // Same-grade differences must not widen.
    static_assert(std::is_same<decltype(v - Vector3<double>()), Vector3<double>>::value,
                  "Vector3 - Vector3 stays a Vector3");
    static_assert(std::is_same<decltype(s - v), Multivector3<double>>::value,
                  "Scalar - Vector3 widens to Multivector3");
}

void test_compound_assignment()
{
    section("compound assignment");

    Vector3<double> v(1, 0, 0);
    v += Vector3<double>(0, 1, 0);
    check_vector(v, 1.0, 1.0, 0.0, "Vector3 +=");
    v -= Vector3<double>(1, 0, 0);
    check_vector(v, 0.0, 1.0, 0.0, "Vector3 -=");
    v *= 3.0;
    check_vector(v, 0.0, 3.0, 0.0, "Vector3 *= raw");
    v /= Scalar<double>(3.0);
    check_vector(v, 0.0, 1.0, 0.0, "Vector3 /= Scalar");

    Bivector3<double> b(1, 2, 3);
    b += Bivector3<double>(1, 1, 1);
    check_bivector(b, 2.0, 3.0, 4.0, "Bivector3 +=");
    b *= Scalar<double>(2.0);
    check_bivector(b, 4.0, 6.0, 8.0, "Bivector3 *= Scalar");

    Trivector3<double> t(4.0);
    t += Trivector3<double>(2.0);
    check_trivector(t, 6.0, "Trivector3 +=");
    t /= 2.0;
    check_trivector(t, 3.0, "Trivector3 /= raw");

    Scalar<double> s(10.0);
    s -= Scalar<double>(4.0);
    check_scalar(s, 6.0, "Scalar -=");
    s *= Scalar<double>(2.0);
    check_scalar(s, 12.0, "Scalar *=");

    Multivector3<double> m = make_mv(1, 1, 1, 1, 1, 1, 1, 1);
    m += make_mv(1, 1, 1, 1, 1, 1, 1, 1);
    check_mv(m, 2, 2, 2, 2, 2, 2, 2, 2, "Multivector3 +=");
    m *= Scalar<double>(0.5);
    check_mv(m, 1, 1, 1, 1, 1, 1, 1, 1, "Multivector3 *= Scalar");

    Rotor3<double> r(Scalar<double>(1.0), Bivector3<double>(2.0, 0.0, 0.0));
    r += Rotor3<double>(Scalar<double>(1.0), Bivector3<double>(2.0, 0.0, 0.0));
    check_scalar(r.scalar, 2.0, "Rotor3 += scalar part");
    check_bivector(r.bivector, 4.0, 0.0, 0.0, "Rotor3 += bivector part");
}

void test_scalar_on_the_left()
{
    section("raw scalar on the left");

    // These used to work only with the number on the right.
    check_vector(3.0 * Vector3<double>(1, 2, 3), 3.0, 6.0, 9.0, "raw * Vector3");
    check_bivector(3.0 * Bivector3<double>(1, 2, 3), 3.0, 6.0, 9.0, "raw * Bivector3");
    check_trivector(3.0 * Trivector3<double>(2.0), 6.0, "raw * Trivector3");

    // Both spellings must agree.
    const Vector3<double> v(1, 2, 3);
    const Vector3<double> left = 2.0 * v;
    const Vector3<double> right = v * 2.0;
    check_vector(left, right.x, right.y, right.z, "2*v == v*2");

    check_mv(2.0 * make_mv(1, 1, 1, 1, 1, 1, 1, 1), 2, 2, 2, 2, 2, 2, 2, 2, "raw * Multivector3");
    const Rotor3<double> r = 2.0 * Rotor3<double>(Scalar<double>(1.0), Bivector3<double>(1, 0, 0));
    check_scalar(r.scalar, 2.0, "raw * Rotor3");
}

void test_grade_projection()
{
    section("grade projection");

    const Multivector3<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8);
    check_scalar(CliffordCore::grade0(m), 1.0, "grade0");
    check_vector(CliffordCore::grade1(m), 2.0, 3.0, 4.0, "grade1");
    check_bivector(CliffordCore::grade2(m), 5.0, 6.0, 7.0, "grade2");
    check_trivector(CliffordCore::grade3(m), 8.0, "grade3");

    const Rotor3<double> r(Scalar<double>(9.0), Bivector3<double>(1, 2, 3));
    check_scalar(CliffordCore::grade0(r), 9.0, "grade0(Rotor3)");
    check_bivector(CliffordCore::grade2(r), 1.0, 2.0, 3.0, "grade2(Rotor3)");

    // Reassembling the grades must give the original back.
    const Multivector3<double> rebuilt = CliffordCore::grade0(m) + CliffordCore::grade1(m)
                                       + CliffordCore::grade2(m) + CliffordCore::grade3(m);
    check_close(mv_difference(rebuilt, m), 0.0, "grades reassemble into the original");
}

void test_rotor_construction()
{
    section("rotor construction");

    check_scalar(CliffordCore::identity_rotor<double>().scalar, 1.0, "identity rotor scalar");
    check_bivector(CliffordCore::identity_rotor<double>().bivector, 0.0, 0.0, 0.0,
                   "identity rotor bivector");

    // Right-hand rule: +90 degrees about +z takes e1 to e2.
    const Rotor3<double> aboutZ =
        CliffordCore::rotor_from_axis_angle(Vector3<double>(0, 0, 1), kPi / 2.0);
    check_vector(CliffordCore::rotate(Vector3<double>(1, 0, 0), aboutZ), 0.0, 1.0, 0.0,
                 "+90 about z takes e1 to e2", 1e-15);
    check_vector(CliffordCore::rotate(Vector3<double>(0, 1, 0), aboutZ), -1.0, 0.0, 0.0,
                 "+90 about z takes e2 to -e1", 1e-15);
    check_vector(CliffordCore::rotate(Vector3<double>(0, 0, 1), aboutZ), 0.0, 0.0, 1.0,
                 "the axis is fixed", 1e-15);

    // +90 about +x takes e2 to e3.
    const Rotor3<double> aboutX =
        CliffordCore::rotor_from_axis_angle(Vector3<double>(1, 0, 0), kPi / 2.0);
    check_vector(CliffordCore::rotate(Vector3<double>(0, 1, 0), aboutX), 0.0, 0.0, 1.0,
                 "+90 about x takes e2 to e3", 1e-15);

    // The axis need not be unit length.
    const Rotor3<double> unnormalised =
        CliffordCore::rotor_from_axis_angle(Vector3<double>(0, 0, 5), kPi / 2.0);
    check_vector(CliffordCore::rotate(Vector3<double>(1, 0, 0), unnormalised), 0.0, 1.0, 0.0,
                 "axis is normalized internally", 1e-15);

    // A zero angle is the identity; a zero axis degrades to the identity too.
    check_vector(CliffordCore::rotate(Vector3<double>(1, 2, 3),
                     CliffordCore::rotor_from_axis_angle(Vector3<double>(0, 0, 1), 0.0)),
                 1.0, 2.0, 3.0, "zero angle is the identity");
    check_scalar(CliffordCore::rotor_from_axis_angle(Vector3<double>(0, 0, 0), 1.0).scalar, 1.0,
                 "zero axis degrades to the identity");

    // rotor_between must actually carry one direction onto the other.
    const Vector3<double> from(1, 2, 3);
    const Vector3<double> to(-2, 1, 0.5);
    const Rotor3<double> between = CliffordCore::rotor_between(from, to);
    const Vector3<double> carried = CliffordCore::rotate(CliffordCore::normalize(from), between);
    const Vector3<double> target = CliffordCore::normalize(to);
    check_vector(carried, target.x, target.y, target.z, "rotor_between carries from onto to", 1e-14);
    check_scalar(CliffordCore::norm(between), 1.0, "rotor_between returns a unit rotor");

    // Identical directions give no rotation.
    check_vector(CliffordCore::rotate(Vector3<double>(1, 0, 0),
                     CliffordCore::rotor_between(Vector3<double>(1, 0, 0), Vector3<double>(2, 0, 0))),
                 1.0, 0.0, 0.0, "rotor_between of parallel vectors is the identity", 1e-14);

    // Opposite directions are the degenerate case: the plane is ambiguous, but
    // the result must still be a genuine half turn rather than NaN.
    const Rotor3<double> flip =
        CliffordCore::rotor_between(Vector3<double>(1, 0, 0), Vector3<double>(-1, 0, 0));
    check_vector(CliffordCore::rotate(Vector3<double>(1, 0, 0), flip), -1.0, 0.0, 0.0,
                 "antiparallel gives a half turn", 1e-14);
    check_scalar(CliffordCore::norm(flip), 1.0, "the half turn is a unit rotor");

    // slerp interpolates: the endpoints are exact and the midpoint is halfway.
    const Rotor3<double> identity = CliffordCore::identity_rotor<double>();
    check_vector(CliffordCore::rotate(Vector3<double>(1, 0, 0),
                     CliffordCore::slerp(identity, aboutZ, 0.0)),
                 1.0, 0.0, 0.0, "slerp at t=0 is the start", 1e-14);
    check_vector(CliffordCore::rotate(Vector3<double>(1, 0, 0),
                     CliffordCore::slerp(identity, aboutZ, 1.0)),
                 0.0, 1.0, 0.0, "slerp at t=1 is the end", 1e-14);
    const double half = std::sqrt(0.5);
    check_vector(CliffordCore::rotate(Vector3<double>(1, 0, 0),
                     CliffordCore::slerp(identity, aboutZ, 0.5)),
                 half, half, 0.0, "slerp at t=0.5 is a 45 degree turn", 1e-14);
}

void test_mixed_products()
{
    section("mixed products");

    const Vector3<double> v(1, 2, 3);
    const Bivector3<double> b(1, 2, 3);
    const Trivector3<double> t(4.0);
    const Multivector3<double> m = make_mv(1, 1, 2, 3, 1, 2, 3, 4);

    // Every mixed product must agree with the general product on promoted operands.
    check_close(mv_difference(v * b, make_mv(0, v.x, v.y, v.z, 0, 0, 0, 0)
                                   * make_mv(0, 0, 0, 0, b.xy, b.xz, b.yz, 0)), 0.0,
                "Vector3 * Bivector3 matches the general product");
    check_close(mv_difference(b * v, make_mv(0, 0, 0, 0, b.xy, b.xz, b.yz, 0)
                                   * make_mv(0, v.x, v.y, v.z, 0, 0, 0, 0)), 0.0,
                "Bivector3 * Vector3 matches the general product");
    check_close(mv_difference(m * v, m * make_mv(0, v.x, v.y, v.z, 0, 0, 0, 0)), 0.0,
                "Multivector3 * Vector3 matches the general product");

    // A bivector squares to minus its squared norm.
    check_close((b * b).scalar.value, -14.0, "b * b scalar part is -|b|^2");
    check_bivector((b * b).bivector, 0.0, 0.0, 0.0, "b * b has no bivector part");

    // Trivector3 * Trivector3 is handled directly on the type rather than here,
    // because the product is always a pure scalar. It must still agree with the
    // general product's scalar part.
    const Multivector3<double> tSquaredGeneral =
        make_mv(0, 0, 0, 0, 0, 0, 0, t.e123) * make_mv(0, 0, 0, 0, 0, 0, 0, t.e123);
    check_close((t * t).value, tSquaredGeneral.scalar.value,
                "Trivector3 * Trivector3 matches the general product");
    check_close(tSquaredGeneral.trivector.e123, 0.0,
                "the general product of two trivectors has no trivector part");

    // Multiplying by the pseudoscalar is the dual, up to the definition.
    check_close(mv_difference(v * t, make_mv(0, v.x, v.y, v.z, 0, 0, 0, 0)
                                   * make_mv(0, 0, 0, 0, 0, 0, 0, t.e123)), 0.0,
                "Vector3 * Trivector3 matches the general product");

    // Wedge with a bivector: a vector lying in the plane wedges to zero, one
    // outside it gives the pseudoscalar.
    check_trivector(Vector3<double>(1, 0, 0) ^ Bivector3<double>(1, 0, 0), 0.0,
                    "e1 ^ e12 is zero, e1 lies in that plane");
    check_trivector(Vector3<double>(0, 1, 0) ^ Bivector3<double>(1, 0, 0), 0.0,
                    "e2 ^ e12 is zero too");
    check_trivector(Vector3<double>(0, 0, 1) ^ Bivector3<double>(1, 0, 0), 1.0,
                    "e3 ^ e12 is the pseudoscalar");
    // The wedge is graded-commutative for grade 1 with grade 2.
    check_trivector(Bivector3<double>(1, 0, 0) ^ Vector3<double>(0, 0, 1), 1.0,
                    "e12 ^ e3 equals e3 ^ e12");
    check_trivector(CliffordCore::wedge_product(v, b), (v ^ b).e123,
                    "named wedge_product matches the operator");

    // The payoff: a rotation can now be written the way the textbooks do.
    const Rotor3<double> r =
        CliffordCore::rotor_from_axis_angle(Vector3<double>(0, 0, 1), 1.0);
    const Multivector3<double> spelled = r * v * CliffordCore::reverse(r);
    const Vector3<double> viaSandwich = CliffordCore::sandwich(v, r);
    check_vector(spelled.vector, viaSandwich.x, viaSandwich.y, viaSandwich.z,
                 "r * v * reverse(r) matches sandwich()", 1e-14);
    check_close(spelled.scalar.value, 0.0, "r * v * reverse(r) has no scalar part", 1e-14);
    check_close(spelled.trivector.e123, 0.0, "r * v * reverse(r) has no trivector part", 1e-14);
    check_bivector(spelled.bivector, 0.0, 0.0, 0.0, "r * v * reverse(r) has no bivector part", 1e-14);
}

} // namespace

int main()
{
    test_scalar();
    test_vector3();
    test_bivector3();
    test_trivector3();
    test_multivector3();
    test_rotor3();

    test_dot_product();
    test_wedge_product();
    test_geometric_product();
    test_norm();
    test_inverse();
    test_reverse();
    test_dual();
    test_exp_log();

    test_basis_multiplication_table();
    test_multivector_product();
    test_norm_all_grades();
    test_reverse_all_grades();
    test_inverse_all_grades();
    test_normalize();
    test_rotation();
    test_mixed_grade_addition();
    test_mixed_grade_subtraction();
    test_compound_assignment();
    test_scalar_on_the_left();
    test_grade_projection();
    test_rotor_construction();
    test_mixed_products();

    std::cout << g_checks << " checks, " << g_failures << " failed.\n";
    if (g_failures != 0) {
        std::cout << "TESTS FAILED\n";
        return 1;
    }
    std::cout << "All CliffordCore tests passed\n";
    return 0;
}
