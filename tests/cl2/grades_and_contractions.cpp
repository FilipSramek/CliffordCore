// Cl(2,0) tests -- grades_and_contractions: grade projection, the left and
// right contractions, and the scalar product.
//
// One category of the Cl2 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

void test_grade_projection()
{
    section("grade projection");

    const Multivector<double> m = make_mv(1.0, 2.0, 3.0, 4.0);
    check_scalar(ga::grade0(m), 1.0, "grade0(Multivector)");
    check_vector(ga::grade1(m), 2.0, 3.0, "grade1(Multivector)");
    check_bivector(ga::grade2(m), 4.0, "grade2(Multivector)");

    const Rotor<double> r(Scalar<double>(5.0), Bivector<double>(6.0));
    check_scalar(ga::grade0(r), 5.0, "grade0(Rotor)");
    check_bivector(ga::grade2(r), 6.0, "grade2(Rotor)");

    // The overload Cl(3,0) lacks: grade0 of a Scalar is the identity, which is
    // what lets scalar_product() accept two bivectors.
    check_scalar(ga::grade0(Scalar<double>(7.0)), 7.0, "grade0(Scalar)");
}

void test_contractions()
{
    section("contractions");

    const Vector<double> u(1.0, 2.0), v(3.0, -1.0);
    const Bivector<double> b(2.0), c(-3.0);

    // Each contraction is the matching grade of the general product.
    const Multivector<double> vb = as_mv(v) * make_mv(0, 0, 0, b.xy);
    const Multivector<double> bv = make_mv(0, 0, 0, b.xy) * as_mv(v);
    const Multivector<double> bc = make_mv(0, 0, 0, b.xy) * make_mv(0, 0, 0, c.xy);

    check_scalar(ga::left_contraction(u, v), (as_mv(u) * as_mv(v)).scalar.value, "left(V, V)");
    check_vector_close(ga::left_contraction(v, b), vb.vector, "left(V, B)");
    check_scalar(ga::left_contraction(b, c), bc.scalar.value, "left(B, B)");
    check_vector_close(v | b, vb.vector, "V | B");
    check_scalar(b | c, bc.scalar.value, "B | B");
    check_scalar(b | c, 6.0, "(2 e12) | (-3 e12) == 6");

    check_scalar(ga::right_contraction(u, v), (u | v).value, "right(V, V) == left(V, V)");
    check_vector_close(ga::right_contraction(b, v), bv.vector, "right(B, V)");
    check_scalar(ga::right_contraction(b, c), bc.scalar.value, "right(B, B)");

    // e12 anticommutes with vectors, so right and left differ by a sign.
    check_vector_close(ga::right_contraction(b, v), -ga::left_contraction(v, b),
                       "right(B, V) == -left(V, B)");

    // e1 contracted into e12 is e2: the quarter turn again.
    check_vector(Vector<double>(1, 0) | Bivector<double>(1), 0.0, 1.0, "e1 | e12 == e2");

    // scalar_product works for every pair whose product has a grade 0 part.
    const Rotor<double> r(Scalar<double>(0.6), Bivector<double>(0.8));
    const Rotor<double> s(Scalar<double>(2.0), Bivector<double>(-1.0));
    const Multivector<double> m = make_mv(1.0, 2.0, 3.0, 4.0);
    const Multivector<double> n = make_mv(-1.0, 0.5, 2.0, -3.0);
    check_scalar(ga::scalar_product(u, v), (u | v).value, "scalar_product(V, V)");
    check_scalar(ga::scalar_product(b, c), (b * c).value, "scalar_product(B, B)");
    check_scalar(ga::scalar_product(r, s), (r * s).scalar.value, "scalar_product(R, R)");
    check_scalar(ga::scalar_product(m, n), (m * n).scalar.value, "scalar_product(M, M)");
}

} // namespace

void run_grades_and_contractions()
{
    test_grade_projection();
    test_contractions();
}

} // namespace cctest
