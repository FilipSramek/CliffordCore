// Cl(3,0) tests -- unary_ops: norm, normalize, reverse, the involutions,
// inverse and dual, for every type.
//
// One category of the Cl3 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

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

} // namespace

void run_unary_ops()
{
    test_norm();
    test_inverse();
    test_reverse();
    test_dual();
    test_norm_all_grades();
    test_reverse_all_grades();
    test_inverse_all_grades();
    test_normalize();
    test_involutions();
}

} // namespace cctest
