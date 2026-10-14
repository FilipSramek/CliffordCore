// Cl(3,0) tests -- products: the dot, wedge and geometric products, the
// general multivector product, and every mixed-grade product checked against
// it.
//
// One category of the Cl3 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

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

} // namespace

void run_products()
{
    test_dot_product();
    test_wedge_product();
    test_geometric_product();
    test_multivector_product();
    test_mixed_products();
}

} // namespace cctest
