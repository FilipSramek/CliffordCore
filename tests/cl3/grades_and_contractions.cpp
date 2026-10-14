// Cl(3,0) tests -- grades_and_contractions: grade projection, the left and
// right contractions, and the scalar product.
//
// One category of the Cl3 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

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

    // Where the product is already a pure Scalar, scalar_product must still
    // work. These did not compile until grade0 gained a Scalar overload.
    const Trivector<double> t1(2.0), t2(3.0);
    check_scalar(ga::scalar_product(t1, t2), -6.0, "scalar_product of two trivectors, e123^2 = -1");
    check_scalar(ga::scalar_product(Scalar<double>(2.0), Scalar<double>(5.0)), 10.0,
                 "scalar_product of two scalars");
    check_scalar(ga::grade0(Scalar<double>(7.0)), 7.0, "grade0(Scalar) is the identity");

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

} // namespace

void run_grades_and_contractions()
{
    test_grade_projection();
    test_contractions();
}

} // namespace cctest
