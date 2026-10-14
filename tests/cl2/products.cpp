// Cl(2,0) tests -- products: the general multivector product, the operator*
// return type of all 25 operand pairs, mixed-grade products checked against
// the general product, and dot and wedge.
//
// One category of the Cl2 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

// Total absolute component difference; 0 means identical.
double mv_difference(const Multivector<double>& a, const Multivector<double>& b)
{
    return std::fabs(a.scalar.value - b.scalar.value)
         + std::fabs(a.vector.x - b.vector.x)
         + std::fabs(a.vector.y - b.vector.y)
         + std::fabs(a.bivector.xy - b.bivector.xy);
}

// Name an operand pair's product type, for the return-type table.
template<typename A, typename B>
using product_t = decltype(std::declval<A>() * std::declval<B>());

// Detects whether a ^ b compiles, so a deliberately missing overload can be
// asserted missing rather than merely not tested.
template<typename A, typename B, typename = void>
struct has_wedge : std::false_type {};

template<typename A, typename B>
struct has_wedge<A, B, decltype(void(std::declval<A>() ^ std::declval<B>()))> : std::true_type {};

void test_multivector_product()
{
    section("multivector product");

    const Multivector<double> a = make_mv(1.0, -2.0, 3.0, 0.5);
    const Multivector<double> b = make_mv(-0.5, 1.0, 2.0, -3.0);
    const Multivector<double> c = make_mv(2.0, 0.5, -1.0, 1.5);
    const Multivector<double> one = basis_blade(0);

    check_mv_close((a * b) * c, a * (b * c), "associativity", 1e-12);
    check_mv_close(a * (b + c), a * b + a * c, "left distributivity", 1e-12);
    check_mv_close((a + b) * c, a * c + b * c, "right distributivity", 1e-12);
    check_mv_close(one * a, a, "1 * a == a");
    check_mv_close(a * one, a, "a * 1 == a");
    check(mv_difference(a * b, b * a) > 1e-6, "the product is not commutative");
    check_mv_close(a * b, ga::geometric_product(a, b), "operator* == geometric_product");

    // Vector * Vector has a direct implementation; it must agree with the table.
    const Vector<double> u(1.0, 2.0), v(3.0, -4.0);
    check_mv_close(u * v, as_mv(u) * as_mv(v), "Vector * Vector == general product");
    check_mv_close(ga::geometric_product(u, v), u * v, "geometric_product(Vector, Vector)");
    const Rotor<double> uv = ga::rotor_product(u, v);
    check_rotor(uv, (u * v).scalar.value, (u * v).bivector.xy, "rotor_product(u, v)");

    // Rotor * Rotor is complex multiplication, also directly implemented.
    const Rotor<double> r(Scalar<double>(0.6), Bivector<double>(-0.8));
    const Rotor<double> s(Scalar<double>(2.0), Bivector<double>(0.5));
    const Multivector<double> rs = ga::to_multivector(r) * ga::to_multivector(s);
    check_rotor(r * s, rs.scalar.value, rs.bivector.xy, "Rotor * Rotor == general product");
    check_close(rs.vector.x, 0.0, "Rotor * Rotor has no vector part (x)");
    check_close(rs.vector.y, 0.0, "Rotor * Rotor has no vector part (y)");
    check_rotor_close(ga::rotor_product(r, s), r * s, "rotor_product(Rotor, Rotor)");

    // The even subalgebra of Cl(2,0) is commutative -- it is the complex
    // numbers. In Cl(3,0) rotors do NOT commute.
    check_rotor_close(r * s, s * r, "rotors commute in two dimensions");
}

void test_operator_return_types()
{
    section("operator* return types");

    using S = Scalar<double>;
    using V = Vector<double>;
    using B = Bivector<double>;
    using M = Multivector<double>;
    using R = Rotor<double>;

    // Every one of the 25 pairs, cell by cell as in docs/cl2.md. The narrow
    // cells are the point of the algebra and must not silently widen.
    static_assert(std::is_same<product_t<S, S>, S>::value, "S*S");
    static_assert(std::is_same<product_t<S, V>, V>::value, "S*V");
    static_assert(std::is_same<product_t<S, B>, B>::value, "S*B");
    static_assert(std::is_same<product_t<S, M>, M>::value, "S*M");
    static_assert(std::is_same<product_t<S, R>, R>::value, "S*R");

    static_assert(std::is_same<product_t<V, S>, V>::value, "V*S");
    static_assert(std::is_same<product_t<V, V>, M>::value, "V*V is a Multivector, as in Cl(3,0)");
    static_assert(std::is_same<product_t<V, B>, V>::value, "V*B is a pure Vector");
    static_assert(std::is_same<product_t<V, M>, M>::value, "V*M");
    static_assert(std::is_same<product_t<V, R>, V>::value, "V*R is a pure Vector");

    static_assert(std::is_same<product_t<B, S>, B>::value, "B*S");
    static_assert(std::is_same<product_t<B, V>, V>::value, "B*V is a pure Vector");
    static_assert(std::is_same<product_t<B, B>, S>::value, "B*B is a pure Scalar");
    static_assert(std::is_same<product_t<B, M>, M>::value, "B*M");
    static_assert(std::is_same<product_t<B, R>, R>::value, "B*R is a Rotor");

    static_assert(std::is_same<product_t<M, S>, M>::value, "M*S");
    static_assert(std::is_same<product_t<M, V>, M>::value, "M*V");
    static_assert(std::is_same<product_t<M, B>, M>::value, "M*B");
    static_assert(std::is_same<product_t<M, M>, M>::value, "M*M");
    static_assert(std::is_same<product_t<M, R>, M>::value, "M*R");

    static_assert(std::is_same<product_t<R, S>, R>::value, "R*S");
    static_assert(std::is_same<product_t<R, V>, V>::value, "R*V is a pure Vector");
    static_assert(std::is_same<product_t<R, B>, R>::value, "R*B is a Rotor");
    static_assert(std::is_same<product_t<R, M>, M>::value, "R*M");
    static_assert(std::is_same<product_t<R, R>, R>::value, "R*R");

    // Vector ^ Bivector would be grade 3, which does not exist here. It is
    // deliberately undefined rather than returning a type that pretends.
    static_assert(has_wedge<V, V>::value, "Vector ^ Vector is defined");
    static_assert(!has_wedge<V, B>::value, "Vector ^ Bivector must NOT be defined");
    static_assert(!has_wedge<B, V>::value, "Bivector ^ Vector must NOT be defined");

    check(true, "all 25 operator* return types and the wedge gap hold");
}

void test_mixed_products()
{
    section("mixed products");

    const Vector<double> v(3.0, -2.0);
    const Bivector<double> b(1.5);
    const Rotor<double> r(Scalar<double>(0.8), Bivector<double>(-0.6));
    const Multivector<double> m = make_mv(0.5, -1.0, 2.0, 0.25);

    const Multivector<double> vb = as_mv(v) * make_mv(0, 0, 0, b.xy);
    const Multivector<double> bv = make_mv(0, 0, 0, b.xy) * as_mv(v);
    const Multivector<double> vr = as_mv(v) * ga::to_multivector(r);
    const Multivector<double> rv = ga::to_multivector(r) * as_mv(v);
    const Multivector<double> br = make_mv(0, 0, 0, b.xy) * ga::to_multivector(r);
    const Multivector<double> rb = ga::to_multivector(r) * make_mv(0, 0, 0, b.xy);

    // Each narrow product equals the general product, AND the grades it drops
    // are zero there -- the narrowing loses nothing.
    check_vector_close(v * b, vb.vector, "Vector * Bivector");
    check_close(vb.scalar.value + vb.bivector.xy, 0.0, "Vector * Bivector drops nothing");
    check_vector_close(b * v, bv.vector, "Bivector * Vector");
    check_close(bv.scalar.value + bv.bivector.xy, 0.0, "Bivector * Vector drops nothing");
    check_vector_close(v * r, vr.vector, "Vector * Rotor");
    check_close(vr.scalar.value + vr.bivector.xy, 0.0, "Vector * Rotor drops nothing");
    check_vector_close(r * v, rv.vector, "Rotor * Vector");
    check_close(rv.scalar.value + rv.bivector.xy, 0.0, "Rotor * Vector drops nothing");
    check_rotor(b * r, br.scalar.value, br.bivector.xy, "Bivector * Rotor");
    check_close(br.vector.x + br.vector.y, 0.0, "Bivector * Rotor drops nothing");
    check_rotor(r * b, rb.scalar.value, rb.bivector.xy, "Rotor * Bivector");
    check_close(rb.vector.x + rb.vector.y, 0.0, "Rotor * Bivector drops nothing");

    // The pseudoscalar anticommutes with vectors but commutes with rotors.
    check_vector_close(b * v, -(v * b), "Bivector * Vector == -(Vector * Bivector)");
    check_rotor_close(b * r, r * b, "Bivector * Rotor == Rotor * Bivector");

    // Everything with a Multivector goes straight through the table.
    check_mv_close(m * v, m * as_mv(v), "Multivector * Vector");
    check_mv_close(v * m, as_mv(v) * m, "Vector * Multivector");
    check_mv_close(m * b, m * make_mv(0, 0, 0, b.xy), "Multivector * Bivector");
    check_mv_close(b * m, make_mv(0, 0, 0, b.xy) * m, "Bivector * Multivector");
    check_mv_close(m * r, m * ga::to_multivector(r), "Multivector * Rotor");
    check_mv_close(r * m, ga::to_multivector(r) * m, "Rotor * Multivector");

    // Bivector * Bivector (a member) agrees with the table's scalar part.
    const Bivector<double> b2(-4.0);
    check_scalar(b * b2, (make_mv(0, 0, 0, b.xy) * make_mv(0, 0, 0, b2.xy)).scalar.value,
                 "Bivector * Bivector == general product");
}

void test_dot_and_wedge()
{
    section("dot and wedge");

    const Vector<double> a(1.0, 2.0), b(3.0, 4.0);
    check_scalar(a | b, 11.0, "(1,2) | (3,4)");
    check_scalar(ga::dot_product(a, b), 11.0, "dot_product");
    check_bivector(a ^ b, -2.0, "(1,2) ^ (3,4)");
    check_bivector(ga::wedge_product(a, b), -2.0, "wedge_product");
    check_bivector(b ^ a, 2.0, "wedge is antisymmetric");
    check_bivector(a ^ a, 0.0, "v ^ v == 0");
    check_bivector(Vector<double>(1, 0) ^ Vector<double>(0, 1), 1.0, "e1 ^ e2 == +e12");

    // The geometric product of two vectors is exactly dot + wedge.
    const Multivector<double> ab = a * b;
    check_mv(ab, 11.0, 0.0, 0.0, -2.0, "a * b == (a|b) + (a^b)");
}

} // namespace

void run_products()
{
    test_multivector_product();
    test_operator_return_types();
    test_mixed_products();
    test_dot_and_wedge();
}

} // namespace cctest
