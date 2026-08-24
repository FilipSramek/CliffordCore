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
void check_scalar(const Scalar<double>& s, double value, const std::string& what)
{
    check_close(s.value, value, what + ".value");
}

void check_vector(const Vector3<double>& v, double x, double y, double z, const std::string& what)
{
    check_close(v.x, x, what + ".x");
    check_close(v.y, y, what + ".y");
    check_close(v.z, z, what + ".z");
}

void check_bivector(const Bivector3<double>& b, double xy, double xz, double yz, const std::string& what)
{
    check_close(b.xy, xy, what + ".xy");
    check_close(b.xz, xz, what + ".xz");
    check_close(b.yz, yz, what + ".yz");
}

void check_trivector(const Trivector3<double>& t, double e123, const std::string& what)
{
    check_close(t.e123, e123, what + ".e123");
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
    check_trivector(a * b, 18.0, "operator* (Trivector3)");
    check_trivector(a / b, 2.0, "operator/ (Trivector3)");
    check_trivector(-a, -6.0, "unary operator-");

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

    // NOTE: dual() is currently a straight component copy, not multiplication by
    // the pseudoscalar. These assertions pin the implemented behavior; if the
    // definition changes to the algebraic dual, update them deliberately.
    const Vector3<double> v(1.0, 2.0, 3.0);
    check_bivector(CliffordCore::dual(v), 1.0, 2.0, 3.0, "dual(Vector3)");

    const Bivector3<double> b(4.0, 5.0, 6.0);
    check_vector(CliffordCore::dual(b), 4.0, 5.0, 6.0, "dual(Bivector3)");

    // Round trip returns the original components.
    const Vector3<double> roundTrip = CliffordCore::dual(CliffordCore::dual(v));
    check_vector(roundTrip, v.x, v.y, v.z, "dual(dual(v)) == v");
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

    std::cout << g_checks << " checks, " << g_failures << " failed.\n";
    if (g_failures != 0) {
        std::cout << "TESTS FAILED\n";
        return 1;
    }
    std::cout << "All CliffordCore tests passed.\n";
    return 0;
}
