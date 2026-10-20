// Cl(2,0) tests -- unary_ops: norm, normalize, the involutions, inverse and
// dual, for every type.
//
// One category of the Cl2 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

// ---------------------------------------------------------------------------
// Magnitudes and involutions
// ---------------------------------------------------------------------------

void test_norms()
{
    section("norms");

    check_scalar(ga::norm(Scalar<double>(-2.0)), 2.0, "norm(Scalar)");
    check_scalar(ga::squared_norm(Scalar<double>(-2.0)), 4.0, "squared_norm(Scalar)");
    check_scalar(ga::norm(Vector<double>(3.0, 4.0)), 5.0, "norm(Vector)");
    check_scalar(ga::squared_norm(Vector<double>(3.0, 4.0)), 25.0, "squared_norm(Vector)");
    check_scalar(ga::norm(Bivector<double>(-6.0)), 6.0, "norm(Bivector)");
    check_scalar(ga::squared_norm(Bivector<double>(-6.0)), 36.0, "squared_norm(Bivector)");
    check_scalar(ga::norm(Rotor<double>(Scalar<double>(3.0), Bivector<double>(4.0))), 5.0, "norm(Rotor)");
    check_scalar(ga::squared_norm(Rotor<double>(Scalar<double>(3.0), Bivector<double>(4.0))), 25.0,
                 "squared_norm(Rotor)");

    // squared_norm is the grade 0 part of m * reverse(m), and the metric is
    // Euclidean, so every component counts.
    const Multivector<double> m = make_mv(1.0, -2.0, 3.0, 0.5);
    check_scalar(ga::squared_norm(m), 1.0 + 4.0 + 9.0 + 0.25, "squared_norm(Multivector)");
    check_scalar(ga::squared_norm(m), (m * ga::reverse(m)).scalar.value, "squared_norm == <m ~m>_0");
    check_scalar(ga::norm(m), std::sqrt(14.25), "norm(Multivector)");
}

void test_normalize()
{
    section("normalize");

    check_scalar(ga::normalize(Scalar<double>(-3.0)), -1.0, "normalize(Scalar) keeps the sign");
    check_vector(ga::normalize(Vector<double>(3.0, 4.0)), 0.6, 0.8, "normalize(Vector)");
    // A bivector has one component, so its unit version is +e12 or -e12.
    check_bivector(ga::normalize(Bivector<double>(-5.0)), -1.0, "normalize(Bivector) keeps orientation");
    check_rotor(ga::normalize(Rotor<double>(Scalar<double>(3.0), Bivector<double>(4.0))), 0.6, 0.8,
                "normalize(Rotor)");
    check_close(ga::norm(ga::normalize(make_mv(1.0, -2.0, 3.0, 0.5))).value, 1.0,
                "normalize(Multivector) is unit");
}

#ifdef NDEBUG
// Release builds only: with the assert compiled out, normalize() hands a zero
// back unchanged. See the same test in tests/cl3/unary_ops.cpp, and audit
// items A1 and A4 in docs/safety-audit.md.
void test_normalize_zero_release()
{
    section("normalize of zero, release build");

    check(ga::normalize(Scalar<double>()) == Scalar<double>(), "normalize(Scalar 0) returns it unchanged");
    check(ga::normalize(Vector<double>()) == Vector<double>(), "normalize(Vector 0) returns it unchanged");
    check(ga::normalize(Bivector<double>()) == Bivector<double>(), "normalize(Bivector 0) returns it unchanged");
    check(ga::normalize(Multivector<double>()) == Multivector<double>(), "normalize(Multivector 0) returns it unchanged");
    const Rotor<double> zero(Scalar<double>(0.0), Bivector<double>());
    check(ga::normalize(zero) == zero, "normalize(Rotor 0) returns it unchanged");
    const Vector<double> tiny(1e-200, 0.0);
    check(ga::normalize(tiny) == tiny, "normalize(1e-200 e1) returns it unchanged, not unit length");
}
#endif

void test_involutions()
{
    section("involutions");

    const Multivector<double> m = make_mv(1.0, 2.0, 3.0, 4.0);

    // The sign table in involutions.hpp, grade by grade.
    check_mv(ga::reverse(m), 1.0, 2.0, 3.0, -4.0, "reverse: + + -");
    check_mv(ga::involute(m), 1.0, -2.0, -3.0, 4.0, "involute: + - +");
    check_mv(ga::conjugate(m), 1.0, -2.0, -3.0, -4.0, "conjugate: + - -");

    // conjugate is reverse composed with involute, in either order.
    check_mv_close(ga::conjugate(m), ga::reverse(ga::involute(m)), "conjugate == reverse . involute");
    check_mv_close(ga::conjugate(m), ga::involute(ga::reverse(m)), "conjugate == involute . reverse");

    // The per-type overloads agree with the multivector one.
    check_scalar(ga::reverse(Scalar<double>(2.0)), 2.0, "reverse(Scalar)");
    check_vector(ga::reverse(Vector<double>(1.0, 2.0)), 1.0, 2.0, "reverse(Vector)");
    check_bivector(ga::reverse(Bivector<double>(3.0)), -3.0, "reverse(Bivector)");
    check_rotor(ga::reverse(Rotor<double>(Scalar<double>(1.0), Bivector<double>(2.0))), 1.0, -2.0,
                "reverse(Rotor)");
    check_scalar(ga::involute(Scalar<double>(2.0)), 2.0, "involute(Scalar)");
    check_vector(ga::involute(Vector<double>(1.0, 2.0)), -1.0, -2.0, "involute(Vector)");
    check_bivector(ga::involute(Bivector<double>(3.0)), 3.0, "involute(Bivector)");
    check_rotor(ga::involute(Rotor<double>(Scalar<double>(1.0), Bivector<double>(2.0))), 1.0, 2.0,
                "involute(Rotor)");
    check_scalar(ga::conjugate(Scalar<double>(2.0)), 2.0, "conjugate(Scalar)");
    check_vector(ga::conjugate(Vector<double>(1.0, 2.0)), -1.0, -2.0, "conjugate(Vector)");
    check_bivector(ga::conjugate(Bivector<double>(3.0)), -3.0, "conjugate(Bivector)");
    check_rotor(ga::conjugate(Rotor<double>(Scalar<double>(1.0), Bivector<double>(2.0))), 1.0, -2.0,
                "conjugate(Rotor) == reverse(Rotor)");

    // The anti-automorphism and automorphism properties, on full multivectors.
    const Multivector<double> n = make_mv(-0.5, 1.5, -2.0, 0.75);
    check_mv_close(ga::reverse(m * n), ga::reverse(n) * ga::reverse(m), "reverse(ab) == ~b ~a");
    check_mv_close(ga::involute(m * n), ga::involute(m) * ga::involute(n),
                   "involute(ab) == involute(a) involute(b)");
    check_mv_close(ga::conjugate(m * n), ga::conjugate(n) * ga::conjugate(m),
                   "conjugate(ab) == conj(b) conj(a)");
}

void test_inverse()
{
    section("inverse");

    const Multivector<double> one = basis_blade(0);

    check_scalar(Scalar<double>(4.0) * ga::inverse(Scalar<double>(4.0)), 1.0, "s * inverse(s)");

    const Vector<double> v(3.0, -4.0);
    check_scalar((v * ga::inverse(v)).scalar, 1.0, "v * inverse(v)");
    check_vector(ga::inverse(Vector<double>(1, 0)), 1.0, 0.0, "inverse(e1) == e1");

    // e12^2 = -1, so the bivector inverse carries a minus sign.
    const Bivector<double> b(2.5);
    check_scalar(b * ga::inverse(b), 1.0, "b * inverse(b)");
    check_bivector(ga::inverse(Bivector<double>(1.0)), -1.0, "inverse(e12) == -e12");

    // A non-unit rotor, so the division by |r|^2 is exercised.
    const Rotor<double> r(Scalar<double>(3.0), Bivector<double>(4.0));
    check_rotor(r * ga::inverse(r), 1.0, 0.0, "r * inverse(r)");
    check_rotor(ga::inverse(r) * r, 1.0, 0.0, "inverse(r) * r");

    // The fact the multivector inverse rests on: m * conjugate(m) is a pure
    // scalar in Cl(2,0), and it is a^2 - x^2 - y^2 + b^2.
    const Multivector<double> m = make_mv(1.0, -2.0, 3.0, 0.5);
    const Multivector<double> collapsed = m * ga::conjugate(m);
    check_mv(collapsed, 1.0 - 4.0 - 9.0 + 0.25, 0.0, 0.0, 0.0, "m * conjugate(m) is a pure scalar");

    check_mv_close(m * ga::inverse(m), one, "m * inverse(m) == 1");
    check_mv_close(ga::inverse(m) * m, one, "inverse(m) * m == 1");

    // That denominator is indefinite. Here it is positive...
    const Multivector<double> p = make_mv(3.0, 1.0, 1.0, 2.0);
    check_mv_close(p * ga::inverse(p), one, "invertible with a positive denominator");
    check_mv_close(ga::inverse(p) * p, one, "invertible with a positive denominator (left)");

    // ...and a non-zero multivector can still be singular. (1 + e1)(1 - e1) = 0.
    const Multivector<double> singular = make_mv(1.0, 1.0, 0.0, 0.0);
    check_mv_close(singular * ga::conjugate(singular), Multivector<double>(),
                   "(1 + e1) is a zero divisor");
    check_mv_close(ga::inverse(singular), Multivector<double>(), "inverse(1 + e1) is zero");
}

void test_dual()
{
    section("dual");

    // The table in dual.hpp: 1 -> e12, e1 -> e2, e2 -> -e1, e12 -> -1.
    check_bivector(ga::dual(Scalar<double>(1.0)), 1.0, "dual(1) == e12");
    check_vector(ga::dual(Vector<double>(1, 0)), 0.0, 1.0, "dual(e1) == e2");
    check_vector(ga::dual(Vector<double>(0, 1)), -1.0, 0.0, "dual(e2) == -e1");
    check_scalar(ga::dual(Bivector<double>(1.0)), -1.0, "dual(e12) == -1");

    // dual(A) is A * e12, checked against the general product.
    const Bivector<double> I(1.0);
    const Vector<double> v(3.0, -2.0);
    const Multivector<double> m = make_mv(1.0, 2.0, 3.0, 4.0);
    check_vector_close(ga::dual(v), v * I, "dual(v) == v * e12");
    check_mv_close(ga::dual(m), m * I, "dual(m) == m * e12");

    // Applying it twice negates, on every grade. Correct, not a defect.
    check_scalar(ga::dual(ga::dual(Scalar<double>(5.0))), -5.0, "dual(dual(s)) == -s");
    check_vector(ga::dual(ga::dual(v)), -3.0, 2.0, "dual(dual(v)) == -v");
    check_bivector(ga::dual(ga::dual(Bivector<double>(5.0))), -5.0, "dual(dual(b)) == -b");
    check_mv_close(ga::dual(ga::dual(m)), -m, "dual(dual(m)) == -m");
}

} // namespace

void run_unary_ops()
{
    test_norms();
    test_normalize();
#ifdef NDEBUG
    test_normalize_zero_release();
#endif
    test_involutions();
    test_inverse();
    test_dual();
}

} // namespace cctest
