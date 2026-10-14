// Cl(3,0,1) tests -- products: the general multivector product, dot and wedge,
// closure of the even subalgebra, and every mixed-grade product checked
// against the promoted table.
//
// One category of the PGA suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

void test_multivector_product()
{
    section("multivector product");

    const Multivector<double> a = make_mv(1, -2, 3, 0.5, -1, 2, -3, 1, 0.25, -0.5, 4, 1, -1, 2, -2, 0.75);
    const Multivector<double> b = make_mv(-0.5, 1, 2, -3, 4, -1, 0.5, 2, 3, -2, 1, 0.5, 1.5, -1, 1, -2);
    const Multivector<double> c = make_mv(2, 0.5, -1, 1, 1, 3, -2, 0.5, -1, 1, -0.5, -2, 1, 0.5, 3, 1);
    const Multivector<double> one = basis_blade(0);

    check_mv_close((a * b) * c, a * (b * c), "associativity", 1e-9);
    check_mv_close(a * (b + c), a * b + a * c, "left distributivity", 1e-9);
    check_mv_close((a + b) * c, a * c + b * c, "right distributivity", 1e-9);
    check_mv_close(one * a, a, "1 * a == a");
    check_mv_close(a * one, a, "a * 1 == a");
    check(mv_difference(a * b, b * a) > 1e-6, "the product is not commutative");

    // operator* and the named function are the same thing.
    check_mv_close(a * b, ga::geometric_product(a, b), "operator* == geometric_product");
}

void test_dot_and_wedge()
{
    section("dot and wedge");

    const Vector<double> a(1, 2, 3, 4);   // e0, e1, e2, e3
    const Vector<double> b(5, 6, 7, 8);

    // The dot product ignores e0: only the normals matter.
    check_scalar(a | b, 2 * 6 + 3 * 7 + 4 * 8, "operator| ignores e0");
    check_scalar(ga::dot_product(a, b), (a | b).value, "dot_product == operator|");
    check_scalar(Vector<double>(1, 0, 0, 0) | Vector<double>(1, 0, 0, 0), 0.0, "e0 | e0 == 0");
    check_scalar(a | b, (b | a).value, "dot product is symmetric");

    // The wedge is antisymmetric and lexicographic.
    check_bivector(a ^ b, 1 * 6 - 2 * 5, 1 * 7 - 3 * 5, 1 * 8 - 4 * 5,
                   2 * 7 - 3 * 6, 2 * 8 - 4 * 6, 3 * 8 - 4 * 7, "operator^ components");
    check_bivector(b ^ a, -(a ^ b).e01, -(a ^ b).e02, -(a ^ b).e03, -(a ^ b).e12, -(a ^ b).e13, -(a ^ b).e23,
                   "wedge is antisymmetric");
    const Bivector<double> self = a ^ a;
    check(self.e01 == 0.0 && self.e02 == 0.0 && self.e03 == 0.0
          && self.e12 == 0.0 && self.e13 == 0.0 && self.e23 == 0.0, "v ^ v == 0 exactly");
    check_bivector(ga::wedge_product(a, b), (a ^ b).e01, (a ^ b).e02, (a ^ b).e03, (a ^ b).e12, (a ^ b).e13, (a ^ b).e23,
                   "wedge_product == operator^");

    // Both agree with the general product.
    const Multivector<double> product = a * b;
    check_scalar(product.scalar, (a | b).value, "scalar part of v*w is the dot");
    check_bivector(product.bivector, (a ^ b).e01, (a ^ b).e02, (a ^ b).e03, (a ^ b).e12, (a ^ b).e13, (a ^ b).e23,
                   "bivector part of v*w is the wedge");
    check_mv_close(ga::geometric_product(a, b), ga::geometric_product(ga::detail::promote(a), ga::detail::promote(b)),
                   "Vector*Vector agrees with the 16x16 table");

    // Geometry: the wedge of two planes is the line where they meet.
    const Vector<double> x0 = test_plane(1, 0, 0, 0);   // x = 0
    const Vector<double> y0 = test_plane(0, 1, 0, 0);   // y = 0
    const Vector<double> x1 = test_plane(1, 0, 0, -1);  // x = 1
    check_bivector(x0 ^ y0, 0, 0, 0, 1, 0, 0, "meet(x=0, y=0) is the z-axis e12");
    check_bivector(x1 ^ y0, 0, -1, 0, 1, 0, 0, "meet(x=1, y=0) is e12 - e02");
    const Vector<double> z0 = test_plane(0, 0, 1, 0);
    const Vector<double> z1 = test_plane(0, 0, 1, -1);
    check_scalar((z0 ^ z1).magnitude(), 0.0, "parallel planes meet in an ideal line");
    check_bivector(z0 ^ z1, 0, 0, 1, 0, 0, 0, "meet(z=0, z=1) is the ideal line +e03");
}

void test_even_closure()
{
    section("even subalgebra closure");

    const Rotor<double> r(Scalar<double>(0.5), 0.25, -0.75, 1.0);
    const Rotor<double> s(Scalar<double>(-1.0), 2.0, 0.5, -0.25);
    const Translator<double> t(Scalar<double>(1.0), 0.5, -1.0, 2.0);
    const Translator<double> u(Scalar<double>(2.0), -0.25, 3.0, 1.0);
    const Motor<double> m(Scalar<double>(0.75), Bivector<double>(1, -2, 0.5, 3, -1, 2), Quadvector<double>(-0.5));
    const Motor<double> n(Scalar<double>(-1.5), Bivector<double>(0.5, 1, -1, -2, 0.25, 1), Quadvector<double>(2.0));

    // The closed forms must agree with the one multiplication table.
    check_mv_close(ga::to_multivector(r * s), ga::to_multivector(r) * ga::to_multivector(s), "Rotor * Rotor matches the table");
    check_mv_close(ga::to_multivector(t * u), ga::to_multivector(t) * ga::to_multivector(u), "Translator * Translator matches the table");
    check_mv_close(ga::to_multivector(m * n), ga::to_multivector(m) * ga::to_multivector(n), "Motor * Motor matches the table");
    check_mv_close(ga::to_multivector(r * t), ga::to_multivector(r) * ga::to_multivector(t), "Rotor * Translator matches the table");
    check_mv_close(ga::to_multivector(t * r), ga::to_multivector(t) * ga::to_multivector(r), "Translator * Rotor matches the table");
    check_mv_close(ga::to_multivector(m * r), ga::to_multivector(m) * ga::to_multivector(r), "Motor * Rotor matches the table");
    check_mv_close(ga::to_multivector(r * m), ga::to_multivector(r) * ga::to_multivector(m), "Rotor * Motor matches the table");
    check_mv_close(ga::to_multivector(m * t), ga::to_multivector(m) * ga::to_multivector(t), "Motor * Translator matches the table");
    check_mv_close(ga::to_multivector(t * m), ga::to_multivector(t) * ga::to_multivector(m), "Translator * Motor matches the table");

    // Return types: closures stay narrow, everything else is a Multivector.
    static_assert(std::is_same<decltype(r * s), Rotor<double>>::value, "Rotor * Rotor is a Rotor");
    static_assert(std::is_same<decltype(t * u), Translator<double>>::value, "Translator * Translator is a Translator");
    static_assert(std::is_same<decltype(m * n), Motor<double>>::value, "Motor * Motor is a Motor");
    static_assert(std::is_same<decltype(r * t), Motor<double>>::value, "Rotor * Translator is a Motor");
    static_assert(std::is_same<decltype(t * r), Motor<double>>::value, "Translator * Rotor is a Motor");
    static_assert(std::is_same<decltype(m * r), Motor<double>>::value, "Motor * Rotor is a Motor");
    static_assert(std::is_same<decltype(t * m), Motor<double>>::value, "Translator * Motor is a Motor");

    // The named spellings are the same products.
    check_mv_close(ga::to_multivector(ga::rotor_product(r, s)), ga::to_multivector(r * s), "rotor_product == operator*");
    check_mv_close(ga::to_multivector(ga::translator_product(t, u)), ga::to_multivector(t * u), "translator_product == operator*");
    check_mv_close(ga::to_multivector(ga::motor_product(m, n)), ga::to_multivector(m * n), "motor_product == operator*");

    // Normalised translators compose by adding their displacements.
    const Translator<double> tx(Scalar<double>(1.0), -0.5, 0, 0);   // +1 along x
    const Translator<double> ty(Scalar<double>(1.0), 0, -1.0, 0);   // +2 along y
    check_translator(tx * ty, 1, -0.5, -1.0, 0, "translators add their ideal parts");
    check_translator(ty * tx, 1, -0.5, -1.0, 0, "translators commute");

    // Two planes multiply to a motor; motor_product packs it.
    const Vector<double> p = test_plane(1, 0, 0, 0), q = test_plane(0, 1, 0, -1);
    check_mv_close(ga::to_multivector(ga::motor_product(p, q)), p * q, "motor_product(plane, plane) == plane * plane");
    // Two points too: grades 0 and 2 only.
    const Trivector<double> P = test_point(1, 2, 3), Q = test_point(-1, 0, 2);
    const Multivector<double> PQ = P * Q;
    static_assert(std::is_same<decltype(P * Q), Multivector<double>>::value, "Trivector * Trivector is a Multivector in PGA");
    check(PQ.vector.magnitude().value == 0.0 && PQ.vector.e0 == 0.0 && PQ.trivector.e123 == 0.0
          && PQ.trivector.e012 == 0.0 && PQ.trivector.e013 == 0.0 && PQ.trivector.e023 == 0.0
          && PQ.quadvector.e0123 == 0.0, "point * point has grades 0 and 2 only");
    check_scalar(PQ.scalar, -1.0, "unit points multiply to scalar -1 plus an ideal line");
    check_mv_close(ga::to_multivector(ga::motor_product(P, Q)), PQ, "motor_product(point, point) == point * point");
}

void test_mixed_products()
{
    section("mixed products");

    const Vector<double> v(1, -2, 0.5, 3);
    const Bivector<double> b(2, -1, 0.5, 1, 3, -2);
    const Trivector<double> t(-1, 2, 0.5, 1.5);
    const Quadvector<double> q(2.5);
    const Multivector<double> m = make_mv(1, -2, 3, 0.5, -1, 2, -3, 1, 0.25, -0.5, 4, 1, -1, 2, -2, 0.75);
    const Rotor<double> r(Scalar<double>(0.5), 0.25, -0.75, 1.0);
    const Translator<double> tr(Scalar<double>(1.0), 0.5, -1.0, 2.0);
    const Motor<double> mo(Scalar<double>(0.75), Bivector<double>(1, -2, 0.5, 3, -1, 2), Quadvector<double>(-0.5));

    const Multivector<double> V = ga::detail::promote(v), B = ga::detail::promote(b), Tt = ga::detail::promote(t);
    const Multivector<double> Q = ga::detail::promote(q), R = ga::detail::promote(r), Tr = ga::detail::promote(tr);
    const Multivector<double> Mo = ga::detail::promote(mo);

    // A representative sweep of the 58 pairs against the table, both orders.
    check_mv_close(v * b, V * B, "Vector * Bivector");
    check_mv_close(b * v, B * V, "Bivector * Vector");
    check_mv_close(v * t, V * Tt, "Vector * Trivector");
    check_mv_close(t * v, Tt * V, "Trivector * Vector");
    check_mv_close(v * q, V * Q, "Vector * Quadvector");
    check_mv_close(q * v, Q * V, "Quadvector * Vector");
    check_mv_close(b * b, B * B, "Bivector * Bivector");
    check_mv_close(b * t, B * Tt, "Bivector * Trivector");
    check_mv_close(t * b, Tt * B, "Trivector * Bivector");
    check_mv_close(b * q, B * Q, "Bivector * Quadvector");
    check_mv_close(q * b, Q * B, "Quadvector * Bivector");
    check_mv_close(t * t, Tt * Tt, "Trivector * Trivector");
    check_mv_close(t * q, Tt * Q, "Trivector * Quadvector");
    check_mv_close(q * t, Q * Tt, "Quadvector * Trivector");
    check_mv_close(m * v, m * V, "Multivector * Vector");
    check_mv_close(v * m, V * m, "Vector * Multivector");
    check_mv_close(m * b, m * B, "Multivector * Bivector");
    check_mv_close(t * m, Tt * m, "Trivector * Multivector");
    check_mv_close(m * q, m * Q, "Multivector * Quadvector");
    check_mv_close(q * m, Q * m, "Quadvector * Multivector");
    check_mv_close(r * v, R * V, "Rotor * Vector");
    check_mv_close(v * r, V * R, "Vector * Rotor");
    check_mv_close(tr * t, Tr * Tt, "Translator * Trivector");
    check_mv_close(t * tr, Tt * Tr, "Trivector * Translator");
    check_mv_close(mo * b, Mo * B, "Motor * Bivector");
    check_mv_close(b * mo, B * Mo, "Bivector * Motor");
    check_mv_close(mo * q, Mo * Q, "Motor * Quadvector");
    check_mv_close(m * mo, m * Mo, "Multivector * Motor");
    check_mv_close(r * m, R * m, "Rotor * Multivector");
    check_mv_close(tr * m, Tr * m, "Translator * Multivector");

    // A line squares to a study number: scalar plus e0123, nothing else.
    const Multivector<double> bb = b * b;
    check(bb.vector.e0 == 0.0 && bb.vector.e1 == 0.0 && bb.bivector.e01 == 0.0 && bb.bivector.e12 == 0.0
          && bb.trivector.e123 == 0.0, "a bivector squared has grades 0 and 4 only");
    check_scalar(bb.scalar, -(1 * 1 + 3 * 3 + 2 * 2), "the scalar part of B*B is minus the Euclidean norm squared");
    check_quadvector(bb.quadvector, 2 * (2 * -2 - (-1) * 3 + 0.5 * 1), "the e0123 part of B*B is 2(e01 e23 - e02 e13 + e03 e12)");

    // Higher wedges via the table.
    check_trivector(v ^ b, (V * B).trivector.e012, (V * B).trivector.e013, (V * B).trivector.e023, (V * B).trivector.e123,
                    "Vector ^ Bivector is the grade 3 part");
    check_trivector(b ^ v, (B * V).trivector.e012, (B * V).trivector.e013, (B * V).trivector.e023, (B * V).trivector.e123,
                    "Bivector ^ Vector is the grade 3 part");
    check_quadvector(v ^ t, (V * Tt).quadvector.e0123, "Vector ^ Trivector is the grade 4 part");
    check_quadvector(t ^ v, (Tt * V).quadvector.e0123, "Trivector ^ Vector is the grade 4 part");
    check_quadvector(b ^ b, (B * B).quadvector.e0123, "Bivector ^ Bivector is the grade 4 part");
    check_quadvector(b ^ Bivector<double>(0, 0, 0, 1, 0, 0), (B * ga::detail::promote(Bivector<double>(0, 0, 0, 1, 0, 0))).quadvector.e0123,
                     "Bivector ^ Bivector against a second line");
    check_trivector(ga::wedge_product(v, b), (v ^ b).e012, (v ^ b).e013, (v ^ b).e023, (v ^ b).e123, "wedge_product(v, b) == v ^ b");
    check_trivector(ga::wedge_product(b, v), (b ^ v).e012, (b ^ v).e013, (b ^ v).e023, (b ^ v).e123, "wedge_product(b, v) == b ^ v");
    check_quadvector(ga::wedge_product(v, t), (v ^ t).e0123, "wedge_product(v, t) == v ^ t");
    check_quadvector(ga::wedge_product(t, v), (t ^ v).e0123, "wedge_product(t, v) == t ^ v");
    check_quadvector(ga::wedge_product(b, b), (b ^ b).e0123, "wedge_product(b, b) == b ^ b");

    // The wedge of a plane with a line is the point where they meet.
    const Vector<double> z1 = test_plane(0, 0, 1, -1);
    const Bivector<double> zAxis(0, 0, 0, 1, 0, 0);
    const Trivector<double> hit = z1 ^ zAxis;
    check_close(hit.e012 / hit.e123, -1.0, "z=1 meets the z-axis at z = 1 (e012 = -z)");
    check_close(hit.e013 / hit.e123, 0.0, "... with y = 0");
    check_close(hit.e023 / hit.e123, 0.0, "... with x = 0");
}

} // namespace

void run_products()
{
    test_multivector_product();
    test_dot_and_wedge();
    test_even_closure();
    test_mixed_products();
}

} // namespace cctest
