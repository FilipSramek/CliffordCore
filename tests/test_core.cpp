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
#include <sstream>
#include <string>

#include "../include/cliffordcore/cl3/scalar.hpp"
#include "../include/cliffordcore/cl3/vector.hpp"
#include "../include/cliffordcore/cl3/bivector.hpp"
#include "../include/cliffordcore/cl3/trivector.hpp"
#include "../include/cliffordcore/cl3/multivector.hpp"
#include "../include/cliffordcore/cl3/rotor.hpp"
#include "../include/cliffordcore/cl3/operations/dot_product.hpp"
#include "../include/cliffordcore/cl3/operations/wedge_product.hpp"
#include "../include/cliffordcore/cl3/operations/geometric_product.hpp"
#include "../include/cliffordcore/cl3/operations/norm.hpp"
#include "../include/cliffordcore/cl3/operations/inverse.hpp"
#include "../include/cliffordcore/cl3/operations/reverse.hpp"
#include "../include/cliffordcore/cl3/operations/dual.hpp"
#include "../include/cliffordcore/cl3/operations/exp.hpp"
#include "../include/cliffordcore/cl3/operations/log.hpp"
#include "../include/cliffordcore/cl3/operations/normalize.hpp"
#include "../include/cliffordcore/cl3/operations/sandwich.hpp"
#include "../include/cliffordcore/cl3/operations/addition.hpp"
#include "../include/cliffordcore/cl3/operations/subtraction.hpp"
#include "../include/cliffordcore/cl3/operations/mixed_products.hpp"
#include "../include/cliffordcore/cl3/operations/rotor_construction.hpp"
#include "../include/cliffordcore/cl3/operations/grade.hpp"
#include "../include/cliffordcore/cl3/operations/involutions.hpp"
#include "../include/cliffordcore/cl3/operations/contraction.hpp"
#include "../include/cliffordcore/cl3/operations/geometry.hpp"
#include "../include/cliffordcore/cl3/operations/comparison.hpp"
#include "../include/cliffordcore/cl3/operations/stream.hpp"

namespace ga = CliffordCore::Cl3;

using ga::Bivector;
using ga::Multivector;
using ga::Rotor;
using ga::Scalar;
using ga::Trivector;
using ga::Vector;

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

void check_vector(const Vector<double>& v, double x, double y, double z, const std::string& what,
                  double tolerance = 1e-12)
{
    check_close(v.x, x, what + ".x", tolerance);
    check_close(v.y, y, what + ".y", tolerance);
    check_close(v.z, z, what + ".z", tolerance);
}

void check_bivector(const Bivector<double>& b, double xy, double xz, double yz, const std::string& what,
                    double tolerance = 1e-12)
{
    check_close(b.xy, xy, what + ".xy", tolerance);
    check_close(b.xz, xz, what + ".xz", tolerance);
    check_close(b.yz, yz, what + ".yz", tolerance);
}

void check_trivector(const Trivector<double>& t, double e123, const std::string& what,
                     double tolerance = 1e-12)
{
    check_close(t.e123, e123, what + ".e123", tolerance);
}

constexpr double kPi = 3.14159265358979323846;

// Build a multivector from raw components, in grade order.
Multivector<double> make_mv(double s, double x, double y, double z,
                             double xy, double xz, double yz, double t)
{
    return Multivector<double>(Scalar<double>(s), Vector<double>(x, y, z),
                                Bivector<double>(xy, xz, yz), Trivector<double>(t));
}

void check_mv(const Multivector<double>& m, double s, double x, double y, double z,
              double xy, double xz, double yz, double t, const std::string& what)
{
    check_close(m.scalar.value, s, what + ".scalar");
    check_vector(m.vector, x, y, z, what);
    check_bivector(m.bivector, xy, xz, yz, what);
    check_close(m.trivector.e123, t, what + ".trivector");
}

// Total absolute component difference; 0 means identical.
double mv_difference(const Multivector<double>& a, const Multivector<double>& b)
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

Multivector<double> rotor_to_mv(const Rotor<double>& r)
{
    return Multivector<double>(r.scalar, Vector<double>(), r.bivector, Trivector<double>());
}

// ---------------------------------------------------------------------------
// Compile-time checks
//
// Everything in the library is constexpr, so these must evaluate during
// compilation. If constexpr-ness regresses, the build fails rather than the run.
// std::sqrt is not portably constexpr, so norm/magnitude stay out of here.
// ---------------------------------------------------------------------------

constexpr Vector<double> kCtA(1.0, 2.0, 3.0);
constexpr Vector<double> kCtB(4.0, 5.0, 6.0);
static_assert((kCtA | kCtB).value == 32.0, "dot product must be constexpr");
static_assert((kCtA ^ kCtB).xy == -3.0, "wedge product must be constexpr");
static_assert(Scalar<double>(2.0).value == 2.0, "Scalar construction must be constexpr");
static_assert(Vector<double>().x == 0.0, "default Vector must be constexpr and zeroed");

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
    check_vector(two * Vector<double>(1.0, 2.0, 3.0), 2.0, 4.0, 6.0, "Scalar * Vector");
    check_bivector(two * Bivector<double>(1.0, 2.0, 3.0), 2.0, 4.0, 6.0, "Scalar * Bivector");
    check_trivector(two * Trivector<double>(5.0), 10.0, "Scalar * Trivector");

    // Composite types: the scalar must distribute over every grade.
    const Multivector<double> mv(
        Scalar<double>(1.0), Vector<double>(2.0, 3.0, 4.0),
        Bivector<double>(5.0, 6.0, 7.0), Trivector<double>(8.0));
    const Multivector<double> scaledMv = two * mv;
    check_scalar(scaledMv.scalar, 2.0, "Scalar * Multivector scalar");
    check_vector(scaledMv.vector, 4.0, 6.0, 8.0, "Scalar * Multivector vector");
    check_bivector(scaledMv.bivector, 10.0, 12.0, 14.0, "Scalar * Multivector bivector");
    check_trivector(scaledMv.trivector, 16.0, "Scalar * Multivector trivector");

    const Rotor<double> rot(Scalar<double>(1.0), Bivector<double>(2.0, 3.0, 4.0));
    const Rotor<double> scaledRot = two * rot;
    check_scalar(scaledRot.scalar, 2.0, "Scalar * Rotor scalar");
    check_bivector(scaledRot.bivector, 4.0, 6.0, 8.0, "Scalar * Rotor bivector");
}

void test_vector3()
{
    section("Vector");

    check_vector(Vector<double>(), 0.0, 0.0, 0.0, "default ctor zeroes");
    check_vector(Vector<double>(1.0, 2.0, 3.0), 1.0, 2.0, 3.0, "component ctor");

    const Vector<double> a(1.0, 2.0, 3.0);
    const Vector<double> b(4.0, 5.0, 6.0);

    check_vector(a + b, 5.0, 7.0, 9.0, "operator+");
    check_vector(a - b, -3.0, -3.0, -3.0, "operator-");
    check_vector(-a, -1.0, -2.0, -3.0, "unary operator-");
    check_vector(a * Scalar<double>(2.0), 2.0, 4.0, 6.0, "operator* (Scalar)");
    check_vector(a / Scalar<double>(2.0), 0.5, 1.0, 1.5, "operator/ (Scalar)");

    // Scalar multiplication agrees across the two spellings.
    const Vector<double> lhs = Scalar<double>(3.0) * a;
    const Vector<double> rhs = a * Scalar<double>(3.0);
    check_vector(lhs, rhs.x, rhs.y, rhs.z, "Scalar*Vector == Vector*Scalar");
}

void test_bivector3()
{
    section("Bivector");

    check_bivector(Bivector<double>(), 0.0, 0.0, 0.0, "default ctor zeroes");
    check_bivector(Bivector<double>(1.0, 2.0, 3.0), 1.0, 2.0, 3.0, "component ctor");

    const Bivector<double> a(1.0, 2.0, 3.0);
    const Bivector<double> b(4.0, 5.0, 6.0);

    check_bivector(a + b, 5.0, 7.0, 9.0, "operator+");
    check_bivector(a - b, -3.0, -3.0, -3.0, "operator-");
    check_bivector(-a, -1.0, -2.0, -3.0, "unary operator-");
    check_bivector(a * Scalar<double>(2.0), 2.0, 4.0, 6.0, "operator* (Scalar)");
    check_bivector(a / Scalar<double>(2.0), 0.5, 1.0, 1.5, "operator/ (Scalar)");

    check_scalar(a.magnitude(), std::sqrt(14.0), "magnitude");
    check_scalar(Bivector<double>().magnitude(), 0.0, "magnitude of zero bivector");

    // magnitude() must agree with the free norm() function.
    check_close(a.magnitude().value, ga::norm(a).value, "magnitude == norm");
}

void test_trivector3()
{
    section("Trivector");

    check_trivector(Trivector<double>(), 0.0, "default ctor zeroes");
    check_trivector(Trivector<double>(7.0), 7.0, "value ctor");

    const Trivector<double> a(6.0);
    const Trivector<double> b(3.0);

    check_trivector(a + b, 9.0, "operator+");
    check_trivector(a - b, 3.0, "operator-");
    check_trivector(-a, -6.0, "unary operator-");

    // The geometric product of two pseudoscalars is a SCALAR, not a trivector:
    // e123 squares to -1, so (6 e123)(3 e123) = -18. Division is a * inverse(b),
    // where the two minus signs cancel and leave a plain 6/3.
    static_assert(std::is_same<decltype(a * b), Scalar<double>>::value,
                  "Trivector * Trivector produces a Scalar");
    static_assert(std::is_same<decltype(a / b), Scalar<double>>::value,
                  "Trivector / Trivector produces a Scalar");
    check_scalar(a * b, -18.0, "operator* (Trivector) is -ab");
    check_scalar(a / b, 2.0, "operator/ (Trivector)");

    // The defining case, and the one the old component-wise version got wrong.
    check_scalar(Trivector<double>(1.0) * Trivector<double>(1.0), -1.0, "e123 * e123 == -1");

    // Division must undo multiplication: (a/b) scaled back by b returns a.
    check_trivector((a / b) * b, a.e123, "(a/b) * b == a");

    check_trivector(a * Scalar<double>(2.0), 12.0, "operator* (Scalar)");
    check_trivector(a / Scalar<double>(2.0), 3.0, "operator/ (Scalar)");

    // Raw numeric overloads exist so `t * 2.0` stays unambiguous now that both
    // Trivector and Scalar are constructible from one T. They must agree with
    // the Scalar spelling.
    check_trivector(a * 2.0, 12.0, "operator* (raw T)");
    check_trivector(a / 2.0, 3.0, "operator/ (raw T)");
    check_close((a * 2.0).e123, (a * Scalar<double>(2.0)).e123, "raw T == Scalar overload");

    check_scalar(a.magnitude(), 6.0, "magnitude");
    check_scalar(Trivector<double>(-6.0).magnitude(), 6.0, "magnitude is non-negative");
}

void test_multivector3()
{
    section("Multivector");

    const Multivector<double> zero;
    check_scalar(zero.scalar, 0.0, "default ctor scalar");
    check_vector(zero.vector, 0.0, 0.0, 0.0, "default ctor vector");
    check_bivector(zero.bivector, 0.0, 0.0, 0.0, "default ctor bivector");
    check_trivector(zero.trivector, 0.0, "default ctor trivector");

    const Multivector<double> a(
        Scalar<double>(1.0), Vector<double>(2.0, 3.0, 4.0),
        Bivector<double>(5.0, 6.0, 7.0), Trivector<double>(8.0));
    const Multivector<double> b(
        Scalar<double>(10.0), Vector<double>(20.0, 30.0, 40.0),
        Bivector<double>(50.0, 60.0, 70.0), Trivector<double>(80.0));

    const Multivector<double> sum = a + b;
    check_scalar(sum.scalar, 11.0, "operator+ scalar");
    check_vector(sum.vector, 22.0, 33.0, 44.0, "operator+ vector");
    check_bivector(sum.bivector, 55.0, 66.0, 77.0, "operator+ bivector");
    check_trivector(sum.trivector, 88.0, "operator+ trivector");

    const Multivector<double> diff = b - a;
    check_scalar(diff.scalar, 9.0, "operator- scalar");
    check_vector(diff.vector, 18.0, 27.0, 36.0, "operator- vector");
    check_bivector(diff.bivector, 45.0, 54.0, 63.0, "operator- bivector");
    check_trivector(diff.trivector, 72.0, "operator- trivector");

    const Multivector<double> negated = -a;
    check_scalar(negated.scalar, -1.0, "unary operator- scalar");
    check_vector(negated.vector, -2.0, -3.0, -4.0, "unary operator- vector");
    check_bivector(negated.bivector, -5.0, -6.0, -7.0, "unary operator- bivector");
    check_trivector(negated.trivector, -8.0, "unary operator- trivector");

    const Multivector<double> scaled = a * Scalar<double>(2.0);
    check_scalar(scaled.scalar, 2.0, "operator* (Scalar) scalar");
    check_vector(scaled.vector, 4.0, 6.0, 8.0, "operator* (Scalar) vector");
    check_bivector(scaled.bivector, 10.0, 12.0, 14.0, "operator* (Scalar) bivector");
    check_trivector(scaled.trivector, 16.0, "operator* (Scalar) trivector");

    const Multivector<double> halved = a / Scalar<double>(2.0);
    check_scalar(halved.scalar, 0.5, "operator/ (Scalar) scalar");
    check_vector(halved.vector, 1.0, 1.5, 2.0, "operator/ (Scalar) vector");
    check_bivector(halved.bivector, 2.5, 3.0, 3.5, "operator/ (Scalar) bivector");
    check_trivector(halved.trivector, 4.0, "operator/ (Scalar) trivector");

    // Scaling agrees whichever side the scalar sits on.
    const Multivector<double> fromLeft = Scalar<double>(2.0) * a;
    check_close(fromLeft.scalar.value, scaled.scalar.value, "Scalar*MV == MV*Scalar scalar");
    check_vector(fromLeft.vector, scaled.vector.x, scaled.vector.y, scaled.vector.z,
                 "Scalar*MV == MV*Scalar vector");
    check_trivector(fromLeft.trivector, scaled.trivector.e123, "Scalar*MV == MV*Scalar trivector");

    // Negation is subtraction from zero.
    const Multivector<double> viaSubtraction = Multivector<double>() - a;
    check_scalar(viaSubtraction.scalar, negated.scalar.value, "-a == 0 - a scalar");
    check_trivector(viaSubtraction.trivector, negated.trivector.e123, "-a == 0 - a trivector");
}

void test_rotor3()
{
    section("Rotor");

    const Rotor<double> zero;
    check_scalar(zero.scalar, 0.0, "default ctor scalar");
    check_bivector(zero.bivector, 0.0, 0.0, 0.0, "default ctor bivector");

    const Rotor<double> a(Scalar<double>(1.0), Bivector<double>(2.0, 3.0, 4.0));
    const Rotor<double> b(Scalar<double>(10.0), Bivector<double>(20.0, 30.0, 40.0));

    const Rotor<double> sum = a + b;
    check_scalar(sum.scalar, 11.0, "operator+ scalar");
    check_bivector(sum.bivector, 22.0, 33.0, 44.0, "operator+ bivector");

    const Rotor<double> diff = b - a;
    check_scalar(diff.scalar, 9.0, "operator- scalar");
    check_bivector(diff.bivector, 18.0, 27.0, 36.0, "operator- bivector");

    const Rotor<double> scaled = a * Scalar<double>(2.0);
    check_scalar(scaled.scalar, 2.0, "operator* (Scalar) scalar");
    check_bivector(scaled.bivector, 4.0, 6.0, 8.0, "operator* (Scalar) bivector");
}

// ---------------------------------------------------------------------------
// Operations
// ---------------------------------------------------------------------------

void test_dot_product()
{
    section("dot_product");

    const Vector<double> a(1.0, 2.0, 3.0);
    const Vector<double> b(4.0, 5.0, 6.0);

    check_scalar(ga::dot_product(a, b), 32.0, "dot_product");
    check_scalar(a | b, 32.0, "operator|");
    check_close((a | b).value, ga::dot_product(a, b).value, "operator| == dot_product");

    // Algebraic properties, not just one hardcoded number.
    check_close((a | b).value, (b | a).value, "dot product is commutative");
    check_scalar(a | a, 14.0, "v|v equals squared norm");

    const Vector<double> x(1.0, 0.0, 0.0);
    const Vector<double> y(0.0, 1.0, 0.0);
    check_scalar(x | y, 0.0, "orthogonal basis vectors have zero dot");
}

void test_wedge_product()
{
    section("wedge_product");

    const Vector<double> a(1.0, 2.0, 3.0);
    const Vector<double> b(4.0, 5.0, 6.0);

    check_bivector(ga::wedge_product(a, b), -3.0, -6.0, -3.0, "wedge_product");
    check_bivector(a ^ b, -3.0, -6.0, -3.0, "operator^");

    // Antisymmetry: a^b == -(b^a).
    const Bivector<double> ab = a ^ b;
    const Bivector<double> ba = b ^ a;
    check_bivector(ba, -ab.xy, -ab.xz, -ab.yz, "wedge is antisymmetric");

    // A vector wedged with itself is zero. Each component is x*y - y*x with
    // identical operands, so this is exactly zero and can be compared with ==.
    const Bivector<double> selfWedge = a ^ a;
    check(selfWedge.xy == 0.0 && selfWedge.xz == 0.0 && selfWedge.yz == 0.0,
          "v^v is exactly zero");

    // Basis check, pinning the (xy, xz, yz) component convention.
    const Vector<double> x(1.0, 0.0, 0.0);
    const Vector<double> y(0.0, 1.0, 0.0);
    check_bivector(x ^ y, 1.0, 0.0, 0.0, "ex^ey is the xy basis bivector");
}

void test_geometric_product()
{
    section("geometric_product");

    const Vector<double> a(1.0, 2.0, 3.0);
    const Vector<double> b(4.0, 5.0, 6.0);

    const Multivector<double> gp = ga::geometric_product(a, b);
    check_scalar(gp.scalar, 32.0, "scalar part is the dot product");
    check_bivector(gp.bivector, -3.0, -6.0, -3.0, "bivector part is the wedge product");
    check_vector(gp.vector, 0.0, 0.0, 0.0, "vector part is zero");
    check_trivector(gp.trivector, 0.0, "trivector part is zero");

    // operator* must agree with the named function.
    const Multivector<double> viaOperator = a * b;
    check_close(viaOperator.scalar.value, gp.scalar.value, "operator* scalar matches");
    check_bivector(viaOperator.bivector, gp.bivector.xy, gp.bivector.xz, gp.bivector.yz,
                   "operator* bivector matches");

    // A vector times itself is a pure scalar equal to its squared norm.
    const Multivector<double> square = a * a;
    check_scalar(square.scalar, 14.0, "v*v scalar part is squared norm");
    check_bivector(square.bivector, 0.0, 0.0, 0.0, "v*v has no bivector part");

    // rotor_product packs the same numbers into a Rotor.
    const Rotor<double> rotor = ga::rotor_product(a, b);
    check_scalar(rotor.scalar, 32.0, "rotor_product scalar");
    check_bivector(rotor.bivector, -3.0, -6.0, -3.0, "rotor_product bivector");
}

void test_norm()
{
    section("norm");

    const Vector<double> a(1.0, 2.0, 3.0);
    check_scalar(ga::norm(a), std::sqrt(14.0), "norm(Vector)");
    check_scalar(ga::squared_norm(a), 14.0, "squared_norm(Vector)");

    const Bivector<double> b(1.0, 2.0, 3.0);
    check_scalar(ga::norm(b), std::sqrt(14.0), "norm(Bivector)");
    check_scalar(ga::squared_norm(b), 14.0, "squared_norm(Bivector)");

    // norm squared must equal squared_norm.
    check_close(ga::norm(a).value * ga::norm(a).value,
                ga::squared_norm(a).value, "norm^2 == squared_norm");

    check_scalar(ga::norm(Vector<double>()), 0.0, "norm of zero vector");
    check_scalar(ga::norm(Vector<double>(3.0, 4.0, 0.0)), 5.0, "3-4-5 triangle");
}

void test_inverse()
{
    section("inverse");

    const Vector<double> a(1.0, 2.0, 3.0);
    const Vector<double> inv = ga::inverse(a);
    check_vector(inv, 1.0 / 14.0, 2.0 / 14.0, 3.0 / 14.0, "inverse(Vector)");

    // Defining property: v * v^-1 == 1, a pure scalar of value one.
    const Multivector<double> product = a * inv;
    check_scalar(product.scalar, 1.0, "v * inverse(v) has unit scalar part");
    check_bivector(product.bivector, 0.0, 0.0, 0.0, "v * inverse(v) has no bivector part");

    // A unit vector is its own inverse.
    const Vector<double> unit(1.0, 0.0, 0.0);
    check_vector(ga::inverse(unit), 1.0, 0.0, 0.0, "unit vector is its own inverse");
}

void test_reverse()
{
    section("reverse");

    const Vector<double> v(1.0, 2.0, 3.0);
    check_vector(ga::reverse(v), 1.0, 2.0, 3.0, "reverse(Vector) is identity");

    const Bivector<double> b(1.0, 2.0, 3.0);
    check_bivector(ga::reverse(b), -1.0, -2.0, -3.0, "reverse(Bivector) negates");

    // Reverse is an involution: applying it twice returns the original.
    const Bivector<double> twice = ga::reverse(ga::reverse(b));
    check_bivector(twice, b.xy, b.xz, b.yz, "reverse(reverse(b)) == b");
}

void test_dual()
{
    section("dual");

    // dual(a) is a * e123. Checked here against the literal basis mapping, and
    // below against the general multivector product.
    const Vector<double> v(1.0, 2.0, 3.0);
    check_bivector(ga::dual(v), 3.0, -2.0, 1.0, "dual(Vector) maps x,y,z to yz,-xz,xy");

    const Bivector<double> b(4.0, 5.0, 6.0);
    check_vector(ga::dual(b), -6.0, 5.0, -4.0, "dual(Bivector)");

    // Grade 0 and grade 3 are dual to each other too.
    check_trivector(ga::dual(Scalar<double>(7.0)), 7.0, "dual(Scalar) is the pseudoscalar");
    check_scalar(ga::dual(Trivector<double>(7.0)), -7.0, "dual(Trivector) negates");

    // Basis checks, pinning each individual mapping.
    check_bivector(ga::dual(Vector<double>(1, 0, 0)), 0.0, 0.0, 1.0, "dual(e1) = e23");
    check_bivector(ga::dual(Vector<double>(0, 1, 0)), 0.0, -1.0, 0.0, "dual(e2) = -e13");
    check_bivector(ga::dual(Vector<double>(0, 0, 1)), 1.0, 0.0, 0.0, "dual(e3) = e12");

    // The pseudoscalar squares to -1, so the dual is an involution up to sign:
    // applying it twice negates rather than returning the original.
    const Vector<double> twice = ga::dual(ga::dual(v));
    check_vector(twice, -v.x, -v.y, -v.z, "dual(dual(v)) == -v");
    check_scalar(ga::dual(ga::dual(Scalar<double>(7.0))), -7.0,
                 "dual(dual(s)) == -s");

    // The definition itself: dual(a) must equal a * e123 under the general
    // product. This is what makes it the algebraic dual rather than a relabelling.
    const Multivector<double> I = make_mv(0, 0, 0, 0, 0, 0, 0, 1);
    const Multivector<double> vAsMv = make_mv(0, v.x, v.y, v.z, 0, 0, 0, 0);
    const Bivector<double> dualV = ga::dual(v);
    check_close(mv_difference(vAsMv * I, make_mv(0, 0, 0, 0, dualV.xy, dualV.xz, dualV.yz, 0)), 0.0,
                "dual(v) == v * e123");

    const Multivector<double> bAsMv = make_mv(0, 0, 0, 0, b.xy, b.xz, b.yz, 0);
    const Vector<double> dualB = ga::dual(b);
    check_close(mv_difference(bAsMv * I, make_mv(0, dualB.x, dualB.y, dualB.z, 0, 0, 0, 0)), 0.0,
                "dual(b) == b * e123");

    // The dual of a wedge is the cross product, which is the classic reason to
    // want a dual in 3D: e1 ^ e2 is the xy plane, whose dual is the e3 axis.
    const Vector<double> e1(1, 0, 0);
    const Vector<double> e2(0, 1, 0);
    check_vector(ga::dual(e1 ^ e2), 0.0, 0.0, -1.0, "dual(e1 ^ e2) is along e3");
}

void test_exp_log()
{
    section("exp/log");

    // exp of the zero bivector is the identity rotor.
    const Rotor<double> identity = ga::exp(Bivector<double>(0.0, 0.0, 0.0));
    check_scalar(identity.scalar, 1.0, "exp(0) scalar");
    check_bivector(identity.bivector, 0.0, 0.0, 0.0, "exp(0) bivector");

    // A quarter turn in the xy plane: cos(pi/2) + sin(pi/2) * xy.
    const Rotor<double> quarter = ga::exp(Bivector<double>(kPi / 2.0, 0.0, 0.0));
    check_close(quarter.scalar.value, 0.0, "exp(pi/2 xy) scalar", 1e-15);
    check_bivector(quarter.bivector, 1.0, 0.0, 0.0, "exp(pi/2 xy) bivector");

    // exp must always produce a unit rotor: scalar^2 + |bivector|^2 == 1.
    const Bivector<double> b(0.3, 0.4, 0.5);
    const Rotor<double> r = ga::exp(b);
    check_close(r.scalar.value * r.scalar.value + ga::squared_norm(r.bivector).value,
                1.0, "exp produces a unit rotor", 1e-15);

    // log of the identity rotor is the zero bivector.
    const Bivector<double> logIdentity =
        ga::log(Rotor<double>(Scalar<double>(1.0), Bivector<double>(0.0, 0.0, 0.0)));
    check_bivector(logIdentity, 0.0, 0.0, 0.0, "log(identity) is zero");

    // log and exp are inverses. Trig costs precision, so allow a looser bound.
    const Bivector<double> there = ga::log(ga::exp(b));
    check_bivector(there, b.xy, b.xz, b.yz, "log(exp(b)) == b");

    const Rotor<double> back = ga::exp(ga::log(r));
    check_close(back.scalar.value, r.scalar.value, "exp(log(r)) scalar", 1e-14);
    check_bivector(back.bivector, r.bivector.xy, r.bivector.xz, r.bivector.yz,
                   "exp(log(r)) bivector");

    // A rotor whose scalar has drifted just outside [-1, 1] must not yield NaN:
    // std::acos would be out of domain without the clamp in log().
    const Bivector<double> drifted =
        ga::log(Rotor<double>(Scalar<double>(1.0 + 2e-16), Bivector<double>(1.0, 0.0, 0.0)));
    check(!std::isnan(drifted.xy) && !std::isnan(drifted.xz) && !std::isnan(drifted.yz),
          "log clamps scalar drift above 1 instead of returning NaN");

    const Bivector<double> driftedLow =
        ga::log(Rotor<double>(Scalar<double>(-1.0 - 2e-16), Bivector<double>(1.0, 0.0, 0.0)));
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
    const Multivector<double> e1  = make_mv(0, 1, 0, 0, 0, 0, 0, 0);
    const Multivector<double> e2  = make_mv(0, 0, 1, 0, 0, 0, 0, 0);
    const Multivector<double> e3  = make_mv(0, 0, 0, 1, 0, 0, 0, 0);
    const Multivector<double> e12 = make_mv(0, 0, 0, 0, 1, 0, 0, 0);
    const Multivector<double> e13 = make_mv(0, 0, 0, 0, 0, 1, 0, 0);
    const Multivector<double> e23 = make_mv(0, 0, 0, 0, 0, 0, 1, 0);
    const Multivector<double> I   = make_mv(0, 0, 0, 0, 0, 0, 0, 1);

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

    const Multivector<double> a = make_mv(1, 2, -3, 4, -5, 6, 7, -8);
    const Multivector<double> b = make_mv(-2, 1, 5, -3, 2, -4, 1, 6);
    const Multivector<double> c = make_mv(3, -1, 2, 1, -2, 3, -1, 2);

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
    const Multivector<double> one = make_mv(1, 0, 0, 0, 0, 0, 0, 0);
    check_close(mv_difference(one * a, a), 0.0, "1 * a == a");
    check_close(mv_difference(a * one, a), 0.0, "a * 1 == a");

    // The general product must reproduce the specialised vector-vector one.
    const Vector<double> u(1, 2, 3);
    const Vector<double> w(4, 5, 6);
    const Multivector<double> viaGeneral =
        make_mv(0, 1, 2, 3, 0, 0, 0, 0) * make_mv(0, 4, 5, 6, 0, 0, 0, 0);
    check_close(mv_difference(viaGeneral, ga::geometric_product(u, w)), 0.0,
                "general product matches geometric_product(v, v)");

    // ...and the specialised rotor-rotor one.
    const Rotor<double> r1(Scalar<double>(0.5), Bivector<double>(0.3, -0.2, 0.7));
    const Rotor<double> r2(Scalar<double>(-0.4), Bivector<double>(0.1, 0.6, -0.5));
    check_close(mv_difference(rotor_to_mv(r1) * rotor_to_mv(r2), rotor_to_mv(r1 * r2)), 0.0,
                "general product matches Rotor * Rotor");

    // A vector times itself is its squared norm, with nothing else left over.
    check_close(mv_difference(make_mv(0, 1, 2, 3, 0, 0, 0, 0) * make_mv(0, 1, 2, 3, 0, 0, 0, 0),
                              make_mv(14, 0, 0, 0, 0, 0, 0, 0)), 0.0, "v * v == |v|^2");
}

void test_norm_all_grades()
{
    section("norm (all grades)");

    check_scalar(ga::norm(Scalar<double>(-3.0)), 3.0, "norm(Scalar) is absolute");
    check_scalar(ga::norm(Trivector<double>(-4.0)), 4.0, "norm(Trivector)");
    check_scalar(ga::squared_norm(Trivector<double>(4.0)), 16.0, "squared_norm(Trivector)");

    // 1 + 4 + 9 + 16 + 25 + 36 + 49 + 64 = 204
    const Multivector<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8);
    check_scalar(ga::squared_norm(m), 204.0, "squared_norm(Multivector)");
    check_scalar(ga::norm(m), std::sqrt(204.0), "norm(Multivector)");

    const Rotor<double> r(Scalar<double>(3.0), Bivector<double>(4.0, 0.0, 0.0));
    check_scalar(ga::squared_norm(r), 25.0, "squared_norm(Rotor)");
    check_scalar(ga::norm(r), 5.0, "norm(Rotor)");

    check_scalar(ga::norm(Multivector<double>()), 0.0, "norm of zero multivector");
}

void test_reverse_all_grades()
{
    section("reverse (all grades)");

    // Reverse flips grade k by (-1)^(k(k-1)/2): grades 0 and 1 keep their sign,
    // grades 2 and 3 are negated.
    check_trivector(ga::reverse(Trivector<double>(5.0)), -5.0, "reverse(Trivector) negates");

    const Multivector<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8);
    check_mv(ga::reverse(m), 1, 2, 3, 4, -5, -6, -7, -8, "reverse(Multivector)");

    // Reverse is an involution on every type.
    check_close(mv_difference(ga::reverse(ga::reverse(m)), m), 0.0,
                "reverse(reverse(m)) == m");
    check_trivector(ga::reverse(ga::reverse(Trivector<double>(5.0))), 5.0,
                    "reverse(reverse(t)) == t");

    // reverse(ab) == reverse(b) reverse(a) -- the defining anti-automorphism.
    const Multivector<double> a = make_mv(1, 2, -3, 4, -5, 6, 7, -8);
    const Multivector<double> b = make_mv(-2, 1, 5, -3, 2, -4, 1, 6);
    check_close(mv_difference(ga::reverse(a * b),
                              ga::reverse(b) * ga::reverse(a)), 0.0,
                "reverse(ab) == reverse(b) reverse(a)", 1e-10);
}

void test_inverse_all_grades()
{
    section("inverse (all grades)");

    check_scalar(ga::inverse(Scalar<double>(4.0)), 0.25, "inverse(Scalar)");

    // Each inverse is checked by its defining property, x * inverse(x) == 1,
    // evaluated with the general product rather than by restating the formula.
    const Multivector<double> one = make_mv(1, 0, 0, 0, 0, 0, 0, 0);

    const Vector<double> v(1.0, 2.0, 3.0);
    check_close(mv_difference(make_mv(0, v.x, v.y, v.z, 0, 0, 0, 0)
                    * make_mv(0, ga::inverse(v).x, ga::inverse(v).y,
                              ga::inverse(v).z, 0, 0, 0, 0), one),
                0.0, "v * inverse(v) == 1");

    // A bivector squares to -|b|^2, so its inverse carries a minus sign.
    const Bivector<double> b(1.0, 2.0, 3.0);
    const Bivector<double> bi = ga::inverse(b);
    check_bivector(bi, -1.0 / 14.0, -2.0 / 14.0, -3.0 / 14.0, "inverse(Bivector) is -b/|b|^2");
    check_close(mv_difference(make_mv(0, 0, 0, 0, b.xy, b.xz, b.yz, 0)
                    * make_mv(0, 0, 0, 0, bi.xy, bi.xz, bi.yz, 0), one),
                0.0, "b * inverse(b) == 1");

    // The pseudoscalar also squares to -1, so the same minus sign applies.
    const Trivector<double> t(4.0);
    const Trivector<double> ti = ga::inverse(t);
    check_trivector(ti, -0.25, "inverse(Trivector) is -t/|t|^2");
    check_close(mv_difference(make_mv(0, 0, 0, 0, 0, 0, 0, t.e123)
                    * make_mv(0, 0, 0, 0, 0, 0, 0, ti.e123), one),
                0.0, "t * inverse(t) == 1");

    // A rotor inverts by reversing, not by dividing itself by its norm.
    const Rotor<double> r(Scalar<double>(0.6), Bivector<double>(0.8, 0.0, 0.0));
    const Rotor<double> ri = ga::inverse(r);
    const Rotor<double> product = r * ri;
    check_scalar(product.scalar, 1.0, "r * inverse(r) scalar is 1");
    check_bivector(product.bivector, 0.0, 0.0, 0.0, "r * inverse(r) has no bivector part");
    // For a unit rotor the inverse is exactly the reverse.
    check_bivector(ri.bivector, ga::reverse(r).bivector.xy,
                   ga::reverse(r).bivector.xz,
                   ga::reverse(r).bivector.yz,
                   "unit rotor inverse == reverse");

    // The general multivector inverse. m/|m|^2 is not it; the conjugate route is.
    const Multivector<double> m = make_mv(1, 2, -3, 4, -5, 6, 7, -8);
    check_close(mv_difference(m * ga::inverse(m), one), 0.0,
                "m * inverse(m) == 1", 1e-10);
    check_close(mv_difference(ga::inverse(m) * m, one), 0.0,
                "inverse(m) * m == 1", 1e-10);

    // Clifford conjugation negates grades 1 and 2 only.
    check_mv(ga::conjugate(make_mv(1, 2, 3, 4, 5, 6, 7, 8)),
             1, -2, -3, -4, -5, -6, -7, 8, "conjugate(m)");
}

void test_normalize()
{
    section("normalize");

    check_vector(ga::normalize(Vector<double>(3.0, 4.0, 0.0)), 0.6, 0.8, 0.0,
                 "normalize(Vector)");
    check_scalar(ga::norm(ga::normalize(Vector<double>(1.0, 2.0, 3.0))), 1.0,
                 "normalized vector has unit norm");
    check_scalar(ga::norm(ga::normalize(Bivector<double>(1.0, 2.0, 3.0))), 1.0,
                 "normalized bivector has unit norm");
    check_trivector(ga::normalize(Trivector<double>(-4.0)), -1.0, "normalize(Trivector)");
    check_scalar(ga::norm(ga::normalize(make_mv(1, 2, 3, 4, 5, 6, 7, 8))), 1.0,
                 "normalized multivector has unit norm");

    // Rotors must be unit length to represent a rotation.
    const Rotor<double> drifted(Scalar<double>(3.0), Bivector<double>(4.0, 0.0, 0.0));
    check_scalar(ga::norm(ga::normalize(drifted)), 1.0,
                 "normalized rotor has unit norm");
    check_scalar(ga::normalize(drifted).scalar, 0.6, "normalize(Rotor) scalar");

    // Normalizing a zero is a precondition violation, not a value: every
    // normalize overload asserts on it. That cannot be checked from inside this
    // harness, since a tripped assert aborts the process, so these cases are
    // deliberately not exercised here.
}

void test_rotation()
{
    section("rotation");

    // exp(theta * B) rotates by 2*theta in the plane B, so pi/4 gives a quarter
    // turn. e12 spans the xy plane, so it turns about z and must leave e3 alone.
    const Rotor<double> quarterXY = ga::exp(Bivector<double>(kPi / 4.0, 0.0, 0.0));
    check_vector(ga::rotate(Vector<double>(1, 0, 0), quarterXY), 0.0, -1.0, 0.0,
                 "e1 rotates to -e2 in the xy plane", 1e-15);
    check_vector(ga::rotate(Vector<double>(0, 0, 1), quarterXY), 0.0, 0.0, 1.0,
                 "e3 is the axis of an xy rotation and is fixed", 1e-15);

    // e23 spans the yz plane, so it turns about x and must leave e1 alone.
    const Rotor<double> quarterYZ = ga::exp(Bivector<double>(0.0, 0.0, kPi / 4.0));
    check_vector(ga::rotate(Vector<double>(1, 0, 0), quarterYZ), 1.0, 0.0, 0.0,
                 "e1 is the axis of a yz rotation and is fixed", 1e-15);

    // Rotation preserves length.
    const Vector<double> v(1.0, 2.0, 3.0);
    check_scalar(ga::norm(ga::rotate(v, quarterXY)), std::sqrt(14.0),
                 "rotation preserves length");

    // The identity rotor leaves everything alone.
    const Rotor<double> identity(Scalar<double>(1.0), Bivector<double>(0.0, 0.0, 0.0));
    check_vector(ga::rotate(v, identity), v.x, v.y, v.z, "identity rotor fixes v");

    // The closed form must agree with the general product R v reverse(R).
    const Rotor<double> r =
        ga::normalize(Rotor<double>(Scalar<double>(0.3), Bivector<double>(0.5, -0.7, 0.2)));
    const Multivector<double> viaGeneral =
        rotor_to_mv(r) * make_mv(0, v.x, v.y, v.z, 0, 0, 0, 0) * rotor_to_mv(ga::reverse(r));
    const Vector<double> viaClosedForm = ga::sandwich(v, r);
    check_vector(viaClosedForm, viaGeneral.vector.x, viaGeneral.vector.y, viaGeneral.vector.z,
                 "sandwich matches R v reverse(R)", 1e-14);
    // The general product must come out a pure vector.
    check_close(viaGeneral.scalar.value, 0.0, "sandwich leaves no scalar part", 1e-14);
    check_close(viaGeneral.trivector.e123, 0.0, "sandwich leaves no trivector part", 1e-14);

    // Rotating twice by R equals rotating once by R*R.
    const Vector<double> twice = ga::rotate(ga::rotate(v, quarterXY), quarterXY);
    const Vector<double> once = ga::rotate(v, quarterXY * quarterXY);
    check_vector(twice, once.x, once.y, once.z, "rotating twice == rotating by the composed rotor", 1e-14);

    // A rotor and its inverse undo each other.
    const Vector<double> roundTrip =
        ga::rotate(ga::rotate(v, quarterXY), ga::inverse(quarterXY));
    check_vector(roundTrip, v.x, v.y, v.z, "rotate then rotate by the inverse returns v", 1e-14);
}

void test_mixed_grade_addition()
{
    section("mixed-grade addition");

    const Scalar<double> s(1.0);
    const Vector<double> v1(1, 0, 0), v2(0, 1, 0), v3(0, 0, 1), v4(1, 1, 1);
    const Bivector<double> b1(1, 0, 0), b2(0, 1, 0), b3(0, 0, 1);
    const Trivector<double> t(5.0);

    // A chain across every grade, in one expression.
    const Multivector<double> m = s + v1 + v2 + v3 + v4 + b1 + b2 + b3 + t;
    check_mv(m, 1, 2, 2, 2, 1, 1, 1, 5, "nine-term mixed chain");

    // Addition is commutative, so the order of the chain must not matter.
    const Multivector<double> reordered = t + b3 + b2 + b1 + v4 + v3 + v2 + v1 + s;
    check_close(mv_difference(m, reordered), 0.0, "chain order does not matter");

    // Every mixed pair lands in the right slot.
    check_mv(s + v1, 1, 1, 0, 0, 0, 0, 0, 0, "Scalar + Vector");
    check_mv(v1 + s, 1, 1, 0, 0, 0, 0, 0, 0, "Vector + Scalar");
    check_mv(s + b1, 1, 0, 0, 0, 1, 0, 0, 0, "Scalar + Bivector");
    check_mv(s + t, 1, 0, 0, 0, 0, 0, 0, 5, "Scalar + Trivector");
    check_mv(v1 + b1, 0, 1, 0, 0, 1, 0, 0, 0, "Vector + Bivector");
    check_mv(v1 + t, 0, 1, 0, 0, 0, 0, 0, 5, "Vector + Trivector");
    check_mv(b1 + t, 0, 0, 0, 0, 1, 0, 0, 5, "Bivector + Trivector");
    check_mv(m + v1, 1, 3, 2, 2, 1, 1, 1, 5, "Multivector + Vector");
    check_mv(v1 + m, 1, 3, 2, 2, 1, 1, 1, 5, "Vector + Multivector");

    // Same-grade sums must not widen.
    static_assert(std::is_same<decltype(v1 + v2), Vector<double>>::value,
                  "Vector + Vector stays a Vector");
    static_assert(std::is_same<decltype(s + b1), Multivector<double>>::value,
                  "Scalar + Bivector widens to Multivector");

    // rotor_sum packs the same numbers as s + b, and to_rotor narrows.
    const Rotor<double> viaSum = ga::rotor_sum(s, b1);
    check_scalar(viaSum.scalar, 1.0, "rotor_sum scalar");
    check_bivector(viaSum.bivector, 1.0, 0.0, 0.0, "rotor_sum bivector");
    const Rotor<double> viaNarrow = ga::to_rotor(s + b1);
    check_scalar(viaNarrow.scalar, 1.0, "to_rotor scalar");
    check_bivector(viaNarrow.bivector, 1.0, 0.0, 0.0, "to_rotor bivector");
    check_mv(ga::to_multivector(viaSum), 1, 0, 0, 0, 1, 0, 0, 0, "to_multivector(rotor)");
}

void test_mixed_grade_subtraction()
{
    section("mixed-grade subtraction");

    const Scalar<double> s(2.0);
    const Vector<double> v(1, 2, 3);
    const Bivector<double> b(4, 5, 6);
    const Trivector<double> t(7.0);
    const Multivector<double> m = make_mv(1, 1, 2, 3, 1, 2, 3, 4);

    check_mv(s - v, 2, -1, -2, -3, 0, 0, 0, 0, "Scalar - Vector");
    check_mv(v - s, -2, 1, 2, 3, 0, 0, 0, 0, "Vector - Scalar");
    check_mv(s - b, 2, 0, 0, 0, -4, -5, -6, 0, "Scalar - Bivector");
    check_mv(s - t, 2, 0, 0, 0, 0, 0, 0, -7, "Scalar - Trivector");
    check_mv(v - b, 0, 1, 2, 3, -4, -5, -6, 0, "Vector - Bivector");
    check_mv(b - v, 0, -1, -2, -3, 4, 5, 6, 0, "Bivector - Vector");
    check_mv(v - t, 0, 1, 2, 3, 0, 0, 0, -7, "Vector - Trivector");
    check_mv(b - t, 0, 0, 0, 0, 4, 5, 6, -7, "Bivector - Trivector");
    check_mv(m - v, 1, 0, 0, 0, 1, 2, 3, 4, "Multivector - Vector");
    check_mv(v - m, -1, 0, 0, 0, -1, -2, -3, -4, "Vector - Multivector");

    // a - b must equal a + (-b), and a - a must vanish.
    check_close(mv_difference(v - b, v + (-b)), 0.0, "a - b == a + (-b)");
    check_close(mv_difference(m - m, Multivector<double>()), 0.0, "m - m == 0");

    // Same-grade differences must not widen.
    static_assert(std::is_same<decltype(v - Vector<double>()), Vector<double>>::value,
                  "Vector - Vector stays a Vector");
    static_assert(std::is_same<decltype(s - v), Multivector<double>>::value,
                  "Scalar - Vector widens to Multivector");
}

void test_compound_assignment()
{
    section("compound assignment");

    Vector<double> v(1, 0, 0);
    v += Vector<double>(0, 1, 0);
    check_vector(v, 1.0, 1.0, 0.0, "Vector +=");
    v -= Vector<double>(1, 0, 0);
    check_vector(v, 0.0, 1.0, 0.0, "Vector -=");
    v *= 3.0;
    check_vector(v, 0.0, 3.0, 0.0, "Vector *= raw");
    v /= Scalar<double>(3.0);
    check_vector(v, 0.0, 1.0, 0.0, "Vector /= Scalar");

    Bivector<double> b(1, 2, 3);
    b += Bivector<double>(1, 1, 1);
    check_bivector(b, 2.0, 3.0, 4.0, "Bivector +=");
    b *= Scalar<double>(2.0);
    check_bivector(b, 4.0, 6.0, 8.0, "Bivector *= Scalar");

    Trivector<double> t(4.0);
    t += Trivector<double>(2.0);
    check_trivector(t, 6.0, "Trivector +=");
    t /= 2.0;
    check_trivector(t, 3.0, "Trivector /= raw");

    Scalar<double> s(10.0);
    s -= Scalar<double>(4.0);
    check_scalar(s, 6.0, "Scalar -=");
    s *= Scalar<double>(2.0);
    check_scalar(s, 12.0, "Scalar *=");

    Multivector<double> m = make_mv(1, 1, 1, 1, 1, 1, 1, 1);
    m += make_mv(1, 1, 1, 1, 1, 1, 1, 1);
    check_mv(m, 2, 2, 2, 2, 2, 2, 2, 2, "Multivector +=");
    m *= Scalar<double>(0.5);
    check_mv(m, 1, 1, 1, 1, 1, 1, 1, 1, "Multivector *= Scalar");

    Rotor<double> r(Scalar<double>(1.0), Bivector<double>(2.0, 0.0, 0.0));
    r += Rotor<double>(Scalar<double>(1.0), Bivector<double>(2.0, 0.0, 0.0));
    check_scalar(r.scalar, 2.0, "Rotor += scalar part");
    check_bivector(r.bivector, 4.0, 0.0, 0.0, "Rotor += bivector part");
}

void test_scalar_on_the_left()
{
    section("raw scalar on the left");

    // These used to work only with the number on the right.
    check_vector(3.0 * Vector<double>(1, 2, 3), 3.0, 6.0, 9.0, "raw * Vector");
    check_bivector(3.0 * Bivector<double>(1, 2, 3), 3.0, 6.0, 9.0, "raw * Bivector");
    check_trivector(3.0 * Trivector<double>(2.0), 6.0, "raw * Trivector");

    // Both spellings must agree.
    const Vector<double> v(1, 2, 3);
    const Vector<double> left = 2.0 * v;
    const Vector<double> right = v * 2.0;
    check_vector(left, right.x, right.y, right.z, "2*v == v*2");

    check_mv(2.0 * make_mv(1, 1, 1, 1, 1, 1, 1, 1), 2, 2, 2, 2, 2, 2, 2, 2, "raw * Multivector");
    const Rotor<double> r = 2.0 * Rotor<double>(Scalar<double>(1.0), Bivector<double>(1, 0, 0));
    check_scalar(r.scalar, 2.0, "raw * Rotor");
}

void test_grade_projection()
{
    section("grade projection");

    const Multivector<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8);
    check_scalar(ga::grade0(m), 1.0, "grade0");
    check_vector(ga::grade1(m), 2.0, 3.0, 4.0, "grade1");
    check_bivector(ga::grade2(m), 5.0, 6.0, 7.0, "grade2");
    check_trivector(ga::grade3(m), 8.0, "grade3");

    const Rotor<double> r(Scalar<double>(9.0), Bivector<double>(1, 2, 3));
    check_scalar(ga::grade0(r), 9.0, "grade0(Rotor)");
    check_bivector(ga::grade2(r), 1.0, 2.0, 3.0, "grade2(Rotor)");

    // Reassembling the grades must give the original back.
    const Multivector<double> rebuilt = ga::grade0(m) + ga::grade1(m)
                                       + ga::grade2(m) + ga::grade3(m);
    check_close(mv_difference(rebuilt, m), 0.0, "grades reassemble into the original");
}

void test_rotor_construction()
{
    section("rotor construction");

    check_scalar(ga::identity_rotor<double>().scalar, 1.0, "identity rotor scalar");
    check_bivector(ga::identity_rotor<double>().bivector, 0.0, 0.0, 0.0,
                   "identity rotor bivector");

    // Right-hand rule: +90 degrees about +z takes e1 to e2.
    const Rotor<double> aboutZ =
        ga::rotor_from_axis_angle(Vector<double>(0, 0, 1), kPi / 2.0);
    check_vector(ga::rotate(Vector<double>(1, 0, 0), aboutZ), 0.0, 1.0, 0.0,
                 "+90 about z takes e1 to e2", 1e-15);
    check_vector(ga::rotate(Vector<double>(0, 1, 0), aboutZ), -1.0, 0.0, 0.0,
                 "+90 about z takes e2 to -e1", 1e-15);
    check_vector(ga::rotate(Vector<double>(0, 0, 1), aboutZ), 0.0, 0.0, 1.0,
                 "the axis is fixed", 1e-15);

    // +90 about +x takes e2 to e3.
    const Rotor<double> aboutX =
        ga::rotor_from_axis_angle(Vector<double>(1, 0, 0), kPi / 2.0);
    check_vector(ga::rotate(Vector<double>(0, 1, 0), aboutX), 0.0, 0.0, 1.0,
                 "+90 about x takes e2 to e3", 1e-15);

    // The axis need not be unit length.
    const Rotor<double> unnormalised =
        ga::rotor_from_axis_angle(Vector<double>(0, 0, 5), kPi / 2.0);
    check_vector(ga::rotate(Vector<double>(1, 0, 0), unnormalised), 0.0, 1.0, 0.0,
                 "axis is normalized internally", 1e-15);

    // A zero angle is the identity; a zero axis degrades to the identity too.
    check_vector(ga::rotate(Vector<double>(1, 2, 3),
                     ga::rotor_from_axis_angle(Vector<double>(0, 0, 1), 0.0)),
                 1.0, 2.0, 3.0, "zero angle is the identity");
    check_scalar(ga::rotor_from_axis_angle(Vector<double>(0, 0, 0), 1.0).scalar, 1.0,
                 "zero axis degrades to the identity");

    // rotor_between must actually carry one direction onto the other.
    const Vector<double> from(1, 2, 3);
    const Vector<double> to(-2, 1, 0.5);
    const Rotor<double> between = ga::rotor_between(from, to);
    const Vector<double> carried = ga::rotate(ga::normalize(from), between);
    const Vector<double> target = ga::normalize(to);
    check_vector(carried, target.x, target.y, target.z, "rotor_between carries from onto to", 1e-14);
    check_scalar(ga::norm(between), 1.0, "rotor_between returns a unit rotor");

    // Identical directions give no rotation.
    check_vector(ga::rotate(Vector<double>(1, 0, 0),
                     ga::rotor_between(Vector<double>(1, 0, 0), Vector<double>(2, 0, 0))),
                 1.0, 0.0, 0.0, "rotor_between of parallel vectors is the identity", 1e-14);

    // Opposite directions are the degenerate case: the plane is ambiguous, but
    // the result must still be a genuine half turn rather than NaN.
    const Rotor<double> flip =
        ga::rotor_between(Vector<double>(1, 0, 0), Vector<double>(-1, 0, 0));
    check_vector(ga::rotate(Vector<double>(1, 0, 0), flip), -1.0, 0.0, 0.0,
                 "antiparallel gives a half turn", 1e-14);
    check_scalar(ga::norm(flip), 1.0, "the half turn is a unit rotor");

    // slerp interpolates: the endpoints are exact and the midpoint is halfway.
    const Rotor<double> identity = ga::identity_rotor<double>();
    check_vector(ga::rotate(Vector<double>(1, 0, 0),
                     ga::slerp(identity, aboutZ, 0.0)),
                 1.0, 0.0, 0.0, "slerp at t=0 is the start", 1e-14);
    check_vector(ga::rotate(Vector<double>(1, 0, 0),
                     ga::slerp(identity, aboutZ, 1.0)),
                 0.0, 1.0, 0.0, "slerp at t=1 is the end", 1e-14);
    const double half = std::sqrt(0.5);
    check_vector(ga::rotate(Vector<double>(1, 0, 0),
                     ga::slerp(identity, aboutZ, 0.5)),
                 half, half, 0.0, "slerp at t=0.5 is a 45 degree turn", 1e-14);
}

void test_mixed_products()
{
    section("mixed products");

    const Vector<double> v(1, 2, 3);
    const Bivector<double> b(1, 2, 3);
    const Trivector<double> t(4.0);
    const Multivector<double> m = make_mv(1, 1, 2, 3, 1, 2, 3, 4);

    // Every mixed product must agree with the general product on promoted operands.
    check_close(mv_difference(v * b, make_mv(0, v.x, v.y, v.z, 0, 0, 0, 0)
                                   * make_mv(0, 0, 0, 0, b.xy, b.xz, b.yz, 0)), 0.0,
                "Vector * Bivector matches the general product");
    check_close(mv_difference(b * v, make_mv(0, 0, 0, 0, b.xy, b.xz, b.yz, 0)
                                   * make_mv(0, v.x, v.y, v.z, 0, 0, 0, 0)), 0.0,
                "Bivector * Vector matches the general product");
    check_close(mv_difference(m * v, m * make_mv(0, v.x, v.y, v.z, 0, 0, 0, 0)), 0.0,
                "Multivector * Vector matches the general product");

    // A bivector squares to minus its squared norm.
    check_close((b * b).scalar.value, -14.0, "b * b scalar part is -|b|^2");
    check_bivector((b * b).bivector, 0.0, 0.0, 0.0, "b * b has no bivector part");

    // Trivector * Trivector is handled directly on the type rather than here,
    // because the product is always a pure scalar. It must still agree with the
    // general product's scalar part.
    const Multivector<double> tSquaredGeneral =
        make_mv(0, 0, 0, 0, 0, 0, 0, t.e123) * make_mv(0, 0, 0, 0, 0, 0, 0, t.e123);
    check_close((t * t).value, tSquaredGeneral.scalar.value,
                "Trivector * Trivector matches the general product");
    check_close(tSquaredGeneral.trivector.e123, 0.0,
                "the general product of two trivectors has no trivector part");

    // Multiplying by the pseudoscalar is the dual, up to the definition.
    check_close(mv_difference(v * t, make_mv(0, v.x, v.y, v.z, 0, 0, 0, 0)
                                   * make_mv(0, 0, 0, 0, 0, 0, 0, t.e123)), 0.0,
                "Vector * Trivector matches the general product");

    // Wedge with a bivector: a vector lying in the plane wedges to zero, one
    // outside it gives the pseudoscalar.
    check_trivector(Vector<double>(1, 0, 0) ^ Bivector<double>(1, 0, 0), 0.0,
                    "e1 ^ e12 is zero, e1 lies in that plane");
    check_trivector(Vector<double>(0, 1, 0) ^ Bivector<double>(1, 0, 0), 0.0,
                    "e2 ^ e12 is zero too");
    check_trivector(Vector<double>(0, 0, 1) ^ Bivector<double>(1, 0, 0), 1.0,
                    "e3 ^ e12 is the pseudoscalar");
    // The wedge is graded-commutative for grade 1 with grade 2.
    check_trivector(Bivector<double>(1, 0, 0) ^ Vector<double>(0, 0, 1), 1.0,
                    "e12 ^ e3 equals e3 ^ e12");
    check_trivector(ga::wedge_product(v, b), (v ^ b).e123,
                    "named wedge_product matches the operator");

    // The payoff: a rotation can now be written the way the textbooks do.
    const Rotor<double> r =
        ga::rotor_from_axis_angle(Vector<double>(0, 0, 1), 1.0);
    const Multivector<double> spelled = r * v * ga::reverse(r);
    const Vector<double> viaSandwich = ga::sandwich(v, r);
    check_vector(spelled.vector, viaSandwich.x, viaSandwich.y, viaSandwich.z,
                 "r * v * reverse(r) matches sandwich()", 1e-14);
    check_close(spelled.scalar.value, 0.0, "r * v * reverse(r) has no scalar part", 1e-14);
    check_close(spelled.trivector.e123, 0.0, "r * v * reverse(r) has no trivector part", 1e-14);
    check_bivector(spelled.bivector, 0.0, 0.0, 0.0, "r * v * reverse(r) has no bivector part", 1e-14);
}


// ---------------------------------------------------------------------------
// The operations added for the initial release
// ---------------------------------------------------------------------------

void test_involutions()
{
    section("involutions");

    const Multivector<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8);

    // Grade involution scales grade k by (-1)^k: odd grades flip.
    check_mv(ga::involute(m), 1, -2, -3, -4, 5, 6, 7, -8, "involute(m)");
    check_scalar(ga::involute(Scalar<double>(3)), 3.0, "involute(Scalar) keeps sign");
    check_vector(ga::involute(Vector<double>(1, 2, 3)), -1.0, -2.0, -3.0,
                 "involute(Vector) flips");
    check_bivector(ga::involute(Bivector<double>(1, 2, 3)), 1.0, 2.0, 3.0,
                   "involute(Bivector) keeps sign");
    check_trivector(ga::involute(Trivector<double>(4)), -4.0, "involute(Trivector) flips");

    // A rotor is purely even, which is what makes it closed under the sandwich.
    const Rotor<double> r(Scalar<double>(1), Bivector<double>(2, 3, 4));
    check_scalar(ga::involute(r).scalar, 1.0, "involute(Rotor) is the identity");
    check_bivector(ga::involute(r).bivector, 2.0, 3.0, 4.0, "involute(Rotor) bivector");

    // An involution: twice is the identity.
    check_close(mv_difference(ga::involute(ga::involute(m)), m), 0.0,
                "involute twice is the identity");

    // The identity that ties the three together, in both compositions.
    check_close(mv_difference(ga::conjugate(m),
                              ga::reverse(ga::involute(m))), 0.0,
                "conjugate == reverse of involute");
    check_close(mv_difference(ga::conjugate(m),
                              ga::involute(ga::reverse(m))), 0.0,
                "conjugate == involute of reverse");

    // conjugate now covers every type, not just Multivector.
    check_vector(ga::conjugate(Vector<double>(1, 2, 3)), -1.0, -2.0, -3.0,
                 "conjugate(Vector)");
    check_trivector(ga::conjugate(Trivector<double>(4)), 4.0,
                    "conjugate(Trivector) keeps sign");
    check_bivector(ga::conjugate(r).bivector, -2.0, -3.0, -4.0, "conjugate(Rotor)");

    // The Scalar overloads the README claimed but that did not exist.
    check_scalar(ga::reverse(Scalar<double>(3)), 3.0, "reverse(Scalar)");
    check_scalar(ga::squared_norm(Scalar<double>(3)), 9.0, "squared_norm(Scalar)");
    check_scalar(ga::normalize(Scalar<double>(-3)), -1.0, "normalize(Scalar) gives the sign");
}

void test_contractions()
{
    section("contractions");

    const Vector<double> a(1, 2, 3);
    const Vector<double> b(4, 5, 6);
    const Bivector<double> B(1, 2, 3);
    const Bivector<double> C(4, 5, 6);
    const Trivector<double> t(2);
    const Trivector<double> u(3);

    // For equal grades the left contraction IS the dot product, so operator|
    // keeps the meaning it always had.
    check_scalar(ga::left_contraction(a, b), 32.0, "V _| V is the dot product");
    check_close((a | b).value, ga::left_contraction(a, b).value, "operator| agrees");

    // Every contraction must equal the correspondingly-graded part of the
    // general product -- that is the definition, so check it rather than
    // restating a formula.
    check_vector(ga::left_contraction(a, B), ga::grade1(a * B).x,
                 ga::grade1(a * B).y, ga::grade1(a * B).z,
                 "V _| B == grade1(a*B)");
    check_bivector(ga::left_contraction(a, t), ga::grade2(a * t).xy,
                   ga::grade2(a * t).xz, ga::grade2(a * t).yz,
                   "V _| T == grade2(a*t)");
    check_close(ga::left_contraction(B, C).value, ga::grade0(B * C).value,
                "B _| B == grade0(B*C)");
    check_close(ga::left_contraction(t, u).value, -6.0, "T _| T is -6, since e123^2 = -1");

    // The defining geometric fact: contracting a vector into a plane it lies in
    // gives the perpendicular direction within that plane.
    const Vector<double> e1(1, 0, 0);
    const Vector<double> e2(0, 1, 0);
    const Vector<double> e3(0, 0, 1);
    check_vector(e1 | (e1 ^ e2), 0.0, 1.0, 0.0, "e1 _| (e1^e2) == e2");
    check_vector(e3 | (e1 ^ e2), 0.0, 0.0, 0.0, "e3 _| (e1^e2) == 0, e3 is outside that plane");

    // Scalar product: the grade 0 part of any product.
    check_scalar(ga::scalar_product(a, b), 32.0, "scalar_product of two vectors");
    check_close(ga::scalar_product(B, C).value, ga::grade0(B * C).value,
                "scalar_product == grade0 of the product");

    // Right contraction is the reverse-mirror of the left one.
    check_vector(ga::right_contraction(B, a),
                 ga::reverse(ga::left_contraction(
                     ga::reverse(a), ga::reverse(B))).x,
                 ga::reverse(ga::left_contraction(
                     ga::reverse(a), ga::reverse(B))).y,
                 ga::reverse(ga::left_contraction(
                     ga::reverse(a), ga::reverse(B))).z,
                 "B |_ a == reverse(reverse(a) _| reverse(B))");
}

void test_geometry()
{
    section("geometry");

    const Vector<double> v(2, 3, 4);
    const Vector<double> n(0, 1, 0);
    const Bivector<double> xy(1, 0, 0);

    // Reflection matches the formula everyone already knows.
    check_vector(ga::reflect(v, n), 2.0, -3.0, 4.0, "reflect in the plane normal to e2");
    const Vector<double> classical = v - 2.0 * (n * (v | n));
    check_vector(ga::reflect(v, n), classical.x, classical.y, classical.z,
                 "reflect matches v - 2(v.n)n");

    // The normal need not be unit length -- reflect divides through by inverse(n).
    check_vector(ga::reflect(v, Vector<double>(0, 5, 0)), 2.0, -3.0, 4.0,
                 "reflect normalizes the plane normal itself");

    // Reflecting twice is the identity.
    check_vector(ga::reflect(ga::reflect(v, n), n), v.x, v.y, v.z,
                 "reflect twice returns the original");

    // A plane can be given as a bivector rather than by its normal.
    check_vector(ga::reflect(v, xy), 2.0, 3.0, -4.0, "reflect in the e12 plane");

    // Projection and rejection split a vector in two, onto a line...
    const Vector<double> u(1, 1, 0);
    const Vector<double> par = ga::project(v, u);
    const Vector<double> perp = ga::reject(v, u);
    check_vector(par, 2.5, 2.5, 0.0, "project onto a line");
    check_vector(par + perp, v.x, v.y, v.z, "project + reject == v, for a line");
    check_scalar(perp | u, 0.0, "the rejection is perpendicular to the line");

    // ...and onto a plane, which is what needed the contraction.
    const Vector<double> inPlane = ga::project(v, xy);
    const Vector<double> outOfPlane = ga::reject(v, xy);
    check_vector(inPlane, 2.0, 3.0, 0.0, "project onto the e12 plane drops the z part");
    check_vector(outOfPlane, 0.0, 0.0, 4.0, "reject from the e12 plane keeps only z");
    check_vector(inPlane + outOfPlane, v.x, v.y, v.z, "project + reject == v, for a plane");

    // The same, on a plane that is not axis aligned.
    const Bivector<double> tilted(1, 2, 3);
    const Vector<double> tiltedSum = ga::project(v, tilted) + ga::reject(v, tilted);
    check_vector(tiltedSum, v.x, v.y, v.z, "project + reject == v, tilted plane", 1e-14);
    const Vector<double> reprojected = ga::project(ga::reject(v, tilted), tilted);
    check_vector(reprojected, 0.0, 0.0, 0.0, "the rejection has nothing left in the plane", 1e-14);
}

void test_sandwich_all_grades()
{
    section("sandwich (all grades)");

    const Rotor<double> r = ga::rotor_from_axis_angle(Vector<double>(0, 0, 1), kPi / 2);

    // Rotating a plane. The e13 plane, turned a quarter turn about z, is e23.
    check_bivector(ga::rotate(Bivector<double>(0, 1, 0), r), 0.0, 0.0, 1.0,
                   "the e13 plane rotates to e23", 1e-15);

    // The rotation plane itself is fixed, exactly as the axis vector is.
    check_bivector(ga::rotate(Bivector<double>(1, 0, 0), r), 1.0, 0.0, 0.0,
                   "the e12 plane is the rotation plane and is fixed", 1e-15);

    // The pseudoscalar commutes with everything, so it cannot be rotated.
    check_trivector(ga::rotate(Trivector<double>(5), r), 5.0,
                    "the pseudoscalar is central and survives untouched", 1e-14);

    // Rotating a multivector must agree grade by grade with rotating the parts.
    const Multivector<double> m(Scalar<double>(1), Vector<double>(1, 0, 0),
                                 Bivector<double>(0, 1, 0), Trivector<double>(5));
    const Multivector<double> rotated = ga::rotate(m, r);
    check_scalar(rotated.scalar, 1.0, "the scalar part is untouched");
    const Vector<double> rv = ga::rotate(Vector<double>(1, 0, 0), r);
    check_vector(rotated.vector, rv.x, rv.y, rv.z, "grade 1 matches rotating the vector alone", 1e-15);
    const Bivector<double> rb = ga::rotate(Bivector<double>(0, 1, 0), r);
    check_bivector(rotated.bivector, rb.xy, rb.xz, rb.yz,
                   "grade 2 matches rotating the bivector alone", 1e-15);

    // The four operator* pairs that were missing, without which none of the
    // above could be written.
    const Bivector<double> b(1, 2, 3);
    const Trivector<double> t(4);
    check_close(mv_difference(r * b, ga::to_multivector(r) * b), 0.0, "Rotor * Bivector");
    check_close(mv_difference(b * r, b * ga::to_multivector(r)), 0.0, "Bivector * Rotor");
    check_close(mv_difference(r * t, ga::to_multivector(r) * t), 0.0, "Rotor * Trivector");
    check_close(mv_difference(t * r, t * ga::to_multivector(r)), 0.0, "Trivector * Rotor");
}

void test_comparison()
{
    section("comparison");

    const Vector<double> a(1, 2, 3);
    const Vector<double> b(1, 2, 3);
    const Vector<double> c(1, 2, 4);

    check(a == b, "identical vectors compare equal");
    check(!(a == c), "differing vectors do not");
    check(a != c, "operator!= is the negation");
    check(!(a != b), "operator!= agrees with operator==");

    check(Scalar<double>(2) == Scalar<double>(2), "Scalar equality");
    check(Bivector<double>(1, 2, 3) == Bivector<double>(1, 2, 3), "Bivector equality");
    check(Trivector<double>(4) == Trivector<double>(4), "Trivector equality");
    check(make_mv(1, 2, 3, 4, 5, 6, 7, 8) == make_mv(1, 2, 3, 4, 5, 6, 7, 8), "Multivector equality");
    check(!(make_mv(1, 2, 3, 4, 5, 6, 7, 8) == make_mv(1, 2, 3, 4, 5, 6, 7, 9)),
          "Multivector inequality catches one differing component");

    // Exact equality is too strict once trigonometry is involved -- which is
    // exactly why approx_equal is a separate function rather than a fuzzy ==.
    const Rotor<double> r = ga::rotor_from_axis_angle(Vector<double>(0, 0, 1), kPi / 2);
    const Vector<double> turned = ga::rotate(Vector<double>(1, 0, 0), r);
    const Vector<double> expected(0, 1, 0);
    check(!(turned == expected), "a rotation leaves residue, so exact == fails");
    check(ga::approx_equal(turned, expected), "approx_equal accepts it");
    check(ga::approx_equal(turned, expected, 1e-9), "approx_equal takes an explicit tolerance");
    check(!ga::approx_equal(turned, expected, 1e-30), "and honours a strict one");

    // A rotor and its negation are the same rotation but different values.
    const Rotor<double> negated(Scalar<double>(-r.scalar.value),
                                 Bivector<double>(-r.bivector.xy, -r.bivector.xz, -r.bivector.yz));
    check(!(r == negated), "a rotor does not equal its negation");
    check(ga::approx_equal(ga::rotate(Vector<double>(1, 0, 0), r),
                                     ga::rotate(Vector<double>(1, 0, 0), negated)),
          "...even though both describe the same rotation");
}

// Exact equality stays usable in a constant expression; a tolerant one could not.
constexpr Vector<double> kCmpA(1, 0, 0);
constexpr Vector<double> kCmpB(1, 0, 0);
static_assert(kCmpA == kCmpB, "operator== is constexpr");
static_assert(!(kCmpA != kCmpB), "operator!= is constexpr");

// ---------------------------------------------------------------------------
// Instantiation sweep
//
// Every check above runs on double and asserts about VALUES. That leaves a
// blind spot: a member of a class template is only compiled when it is called,
// so an entry point nothing calls can be syntactically broken and still ship
// green. That is exactly how to_string() reached main() broken on five of the
// six types while this suite reported "779 checks, 0 failed".
//
// This section closes it. It asserts little -- its job is to CALL every public
// entry point, for float, double and long double, so a compile failure
// anywhere in the API stops the build.
// ---------------------------------------------------------------------------

template <typename T>
void instantiate_every_entry_point()
{
    ga::Scalar<T> s(2);
    ga::Vector<T> v(1, 2, 3);
    ga::Bivector<T> b(1, 2, 3);
    ga::Trivector<T> t(4);
    ga::Multivector<T> m(s, v, b, t);
    ga::Rotor<T> r(ga::Scalar<T>(1), ga::Bivector<T>(0, 0, 0));

    // Same-type arithmetic and unary minus.
    (void)(s + s); (void)(s - s); (void)(s * s); (void)(s / s); (void)(-s);
    (void)(v + v); (void)(v - v); (void)(-v);
    (void)(b + b); (void)(b - b); (void)(-b);
    (void)(t + t); (void)(t - t); (void)(-t); (void)(t * t); (void)(t / t);
    (void)(m + m); (void)(m - m); (void)(-m);
    (void)(r + r); (void)(r - r);

    // Scalar mixing, both spellings and both sides.
    (void)(s * v); (void)(s * b); (void)(s * t); (void)(s * m); (void)(s * r);
    (void)(v * s); (void)(b * s); (void)(t * s); (void)(m * s); (void)(r * s);
    (void)(v / s); (void)(b / s); (void)(t / s); (void)(m / s); (void)(r / s);
    (void)(t * T(2)); (void)(t / T(2));
    (void)(T(2) * v); (void)(T(2) * b); (void)(T(2) * t); (void)(T(2) * m); (void)(T(2) * r);

    // Compound assignment, every operator on every type.
    { ga::Scalar<T> q = s; q += s; q -= s; q *= s; q /= s; }
    { ga::Vector<T> q = v; q += v; q -= v; q *= T(2); q /= T(2); q *= s; q /= s; }
    { ga::Bivector<T> q = b; q += b; q -= b; q *= T(2); q /= T(2); q *= s; q /= s; }
    { ga::Trivector<T> q = t; q += t; q -= t; q *= T(2); q /= T(2); q *= s; q /= s; }
    { ga::Multivector<T> q = m; q += m; q -= m; q *= s; q /= s; }
    { ga::Rotor<T> q = r; q += r; q -= r; q *= s; q /= s; }

    // Products.
    (void)(v | v); (void)ga::dot_product(v, v);
    (void)(v ^ v); (void)ga::wedge_product(v, v);
    (void)(v ^ b); (void)(b ^ v); (void)ga::wedge_product(v, b);
    (void)(v * v); (void)(m * m); (void)(r * r);
    (void)ga::geometric_product(v, v);
    (void)ga::geometric_product(m, m);
    (void)ga::rotor_product(v, v);
    (void)ga::rotor_product(r, r);

    // Every mixed_products.hpp overload, in both orders.
    (void)(v * b); (void)(b * v);
    (void)(b * b);
    (void)(v * t); (void)(t * v);
    (void)(b * t); (void)(t * b);
    (void)(r * v); (void)(v * r);
    (void)(m * v); (void)(v * m);
    (void)(m * b); (void)(b * m);
    (void)(m * t); (void)(t * m);
    (void)(m * r); (void)(r * m);

    // Mixed-grade sums and differences, both orders.
    (void)(s + v); (void)(v + s); (void)(s + b); (void)(b + s); (void)(s + t); (void)(t + s);
    (void)(v + b); (void)(b + v); (void)(v + t); (void)(t + v); (void)(b + t); (void)(t + b);
    (void)(m + s); (void)(s + m); (void)(m + v); (void)(v + m);
    (void)(m + b); (void)(b + m); (void)(m + t); (void)(t + m);
    (void)(s - v); (void)(v - s); (void)(s - b); (void)(b - s); (void)(s - t); (void)(t - s);
    (void)(v - b); (void)(b - v); (void)(v - t); (void)(t - v); (void)(b - t); (void)(t - b);
    (void)(m - s); (void)(s - m); (void)(m - v); (void)(v - m);
    (void)(m - b); (void)(b - m); (void)(m - t); (void)(t - m);

    // Magnitudes, involutions, structure.
    (void)ga::norm(s); (void)ga::norm(v); (void)ga::norm(b);
    (void)ga::norm(t); (void)ga::norm(m); (void)ga::norm(r);
    (void)ga::squared_norm(v); (void)ga::squared_norm(b);
    (void)ga::squared_norm(t); (void)ga::squared_norm(m);
    (void)ga::squared_norm(r);
    (void)ga::normalize(v); (void)ga::normalize(b);
    (void)ga::normalize(t); (void)ga::normalize(m);
    (void)ga::normalize(r);
    (void)ga::reverse(v); (void)ga::reverse(b);
    (void)ga::reverse(t); (void)ga::reverse(m);
    (void)ga::reverse(r);
    (void)ga::inverse(s); (void)ga::inverse(v);
    (void)ga::inverse(b); (void)ga::inverse(t);
    (void)ga::inverse(m); (void)ga::inverse(r);
    (void)ga::conjugate(m);
    (void)ga::dual(s); (void)ga::dual(v);
    (void)ga::dual(b); (void)ga::dual(t);
    (void)ga::grade0(m); (void)ga::grade1(m);
    (void)ga::grade2(m); (void)ga::grade3(m);
    (void)ga::grade0(r); (void)ga::grade2(r);
    (void)ga::to_rotor(m); (void)ga::to_multivector(r);
    (void)ga::rotor_sum(s, b);
    (void)v.magnitude(); (void)b.magnitude(); (void)t.magnitude();

    // Rotations.
    (void)ga::exp(b);
    (void)ga::log(r);
    (void)ga::sandwich(v, r); (void)ga::rotate(v, r);
    (void)ga::identity_rotor<T>();
    (void)ga::rotor_from_axis_angle(v, T(1));
    (void)ga::rotor_between(v, v);
    (void)ga::slerp(r, r, T(0.5));

    // Printing -- the entry points that shipped broken because nothing called them.
    (void)s.to_string(); (void)v.to_string(); (void)b.to_string();
    (void)t.to_string(); (void)m.to_string(); (void)r.to_string();

    // The four operator* pairs added for the bivector sandwich.
    (void)(r * b); (void)(b * r); (void)(r * t); (void)(t * r);

    // Involutions, on every type.
    (void)ga::involute(s); (void)ga::involute(v);
    (void)ga::involute(b); (void)ga::involute(t);
    (void)ga::involute(m); (void)ga::involute(r);
    (void)ga::conjugate(s); (void)ga::conjugate(v);
    (void)ga::conjugate(b); (void)ga::conjugate(t);
    (void)ga::conjugate(m); (void)ga::conjugate(r);
    (void)ga::reverse(s);
    (void)ga::squared_norm(s); (void)ga::normalize(s);

    // Contractions and the scalar product.
    (void)ga::left_contraction(v, v); (void)ga::left_contraction(v, b);
    (void)ga::left_contraction(v, t); (void)ga::left_contraction(b, b);
    (void)ga::left_contraction(b, t); (void)ga::left_contraction(t, t);
    (void)ga::right_contraction(v, v); (void)ga::right_contraction(b, v);
    (void)ga::right_contraction(t, v); (void)ga::right_contraction(b, b);
    (void)ga::right_contraction(t, b); (void)ga::right_contraction(t, t);
    (void)(v | b); (void)(v | t); (void)(b | b); (void)(b | t);
    (void)ga::scalar_product(v, v); (void)ga::scalar_product(b, b);

    // Geometry.
    (void)ga::reflect(v, v); (void)ga::reflect(v, b);
    (void)ga::project(v, v); (void)ga::project(v, b);
    (void)ga::reject(v, v); (void)ga::reject(v, b);

    // Sandwich, on every grade it accepts.
    (void)ga::sandwich(b, r); (void)ga::sandwich(t, r);
    (void)ga::sandwich(m, r);
    (void)ga::rotate(b, r); (void)ga::rotate(t, r);
    (void)ga::rotate(m, r);

    // Comparison, exact and tolerant.
    (void)(s == s); (void)(v == v); (void)(b == b);
    (void)(t == t); (void)(m == m); (void)(r == r);
    (void)(s != s); (void)(v != v); (void)(b != b);
    (void)(t != t); (void)(m != m); (void)(r != r);
    (void)ga::approx_equal(s, s); (void)ga::approx_equal(v, v);
    (void)ga::approx_equal(b, b); (void)ga::approx_equal(t, t);
    (void)ga::approx_equal(m, m); (void)ga::approx_equal(r, r);

    // Stream insertion, for all six.
    {
        std::ostringstream os;
        os << s << v << b << t << m << r;
        (void)os;
    }
}

void test_instantiation_sweep()
{
    section("instantiation sweep");

    instantiate_every_entry_point<float>();
    instantiate_every_entry_point<double>();
    instantiate_every_entry_point<long double>();
    check(true, "every public entry point instantiates for float, double and long double");

    const Scalar<double> s(1.5);
    const Vector<double> v(1, 2, 3);
    const Bivector<double> b(1, 2, 3);
    const Trivector<double> t(4);
    const Multivector<double> m(s, v, b, t);
    const Rotor<double> r(s, b);

    check(!s.to_string().empty(), "Scalar::to_string is non-empty");
    check(!v.to_string().empty(), "Vector::to_string is non-empty");
    check(!b.to_string().empty(), "Bivector::to_string is non-empty");
    check(!t.to_string().empty(), "Trivector::to_string is non-empty");
    check(!m.to_string().empty(), "Multivector::to_string is non-empty");
    check(!r.to_string().empty(), "Rotor::to_string is non-empty");

    // The basis labels must match the component names in docs/conventions.md,
    // or the output misleads about which grade is which.
    check(v.to_string().find("e1") != std::string::npos, "Vector labels its e1 component");
    check(b.to_string().find("e12") != std::string::npos, "Bivector labels its e12 component");
    check(t.to_string().find("e123") != std::string::npos, "Trivector labels its e123 component");

    // Round-trip precision. std::to_string renders this as "0.000000", which is
    // useless for the residuals this library actually produces.
    check(Scalar<double>(1e-17).to_string() != "0.000000",
          "to_string keeps small values instead of flattening them to zero");
    check(Scalar<double>(1e-17).to_string().find("e-17") != std::string::npos,
          "to_string uses scientific notation where it is needed");
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
    test_involutions();
    test_contractions();
    test_geometry();
    test_sandwich_all_grades();
    test_comparison();

    test_instantiation_sweep();

    std::cout << g_checks << " checks, " << g_failures << " failed.\n";
    if (g_failures != 0) {
        std::cout << "TESTS FAILED\n";
        return 1;
    }
    std::cout << "All CliffordCore tests passed\n";
    return 0;
}
