// Cl(3,0,1) tests -- grades_and_contractions: grade projection, the left and
// right contractions, and the scalar product.
//
// One category of the PGA suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

void test_grade_projection()
{
    section("grade projection");

    const Multivector<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
    check_scalar(ga::grade0(m), 1.0, "grade0");
    check_vector(ga::grade1(m), 2, 3, 4, 5, "grade1");
    check_bivector(ga::grade2(m), 6, 7, 8, 9, 10, 11, "grade2");
    check_trivector(ga::grade3(m), 12, 13, 14, 15, "grade3");
    check_quadvector(ga::grade4(m), 16.0, "grade4");
    check_mv_close(ga::grade0(m) + ga::grade1(m) + ga::grade2(m) + ga::grade3(m) + ga::grade4(m), m,
                   "the grades reassemble the multivector");

    const Rotor<double> r(Scalar<double>(1.0), 2, 3, 4);
    check_scalar(ga::grade0(r), 1.0, "grade0(Rotor)");
    check_bivector(ga::grade2(r), 0, 0, 0, 2, 3, 4, "grade2(Rotor) has zero ideal components");

    const Translator<double> t(Scalar<double>(1.0), 2, 3, 4);
    check_scalar(ga::grade0(t), 1.0, "grade0(Translator)");
    check_bivector(ga::grade2(t), 2, 3, 4, 0, 0, 0, "grade2(Translator) has zero Euclidean components");

    const Motor<double> mo(Scalar<double>(1.0), Bivector<double>(2, 3, 4, 5, 6, 7), Quadvector<double>(8.0));
    check_scalar(ga::grade0(mo), 1.0, "grade0(Motor)");
    check_bivector(ga::grade2(mo), 2, 3, 4, 5, 6, 7, "grade2(Motor)");
    check_quadvector(ga::grade4(mo), 8.0, "grade4(Motor)");
}

void test_contractions()
{
    section("contractions");

    const Vector<double> v(1, -2, 0.5, 3);
    const Bivector<double> b(2, -1, 0.5, 1, 3, -2);
    const Trivector<double> t(-1, 2, 0.5, 1.5);
    const Quadvector<double> q(2.5);
    const Vector<double> w(0.5, 1, -1, 2);
    const Bivector<double> c(-1, 0.5, 2, -1, 0.5, 1);
    const Trivector<double> u(2, -1, 1, 0.5);

    // Left contraction is the grade (s - r) part of the product.
    check_scalar(ga::left_contraction(v, w), (v * w).scalar.value, "v _| w is the scalar part");
    check_vector(ga::left_contraction(v, b), (v * b).vector.e0, (v * b).vector.e1, (v * b).vector.e2, (v * b).vector.e3, "v _| B is the vector part");
    check_bivector(ga::left_contraction(v, t), (v * t).bivector.e01, (v * t).bivector.e02, (v * t).bivector.e03,
                   (v * t).bivector.e12, (v * t).bivector.e13, (v * t).bivector.e23, "v _| T is the bivector part");
    check_trivector(ga::left_contraction(v, q), (v * q).trivector.e012, (v * q).trivector.e013, (v * q).trivector.e023,
                    (v * q).trivector.e123, "v _| I is the trivector part");
    check_scalar(ga::left_contraction(b, c), (b * c).scalar.value, "B _| C is the scalar part");
    check_vector(ga::left_contraction(b, t), (b * t).vector.e0, (b * t).vector.e1, (b * t).vector.e2, (b * t).vector.e3, "B _| T is the vector part");
    check_bivector(ga::left_contraction(b, q), (b * q).bivector.e01, (b * q).bivector.e02, (b * q).bivector.e03,
                   (b * q).bivector.e12, (b * q).bivector.e13, (b * q).bivector.e23, "B _| I is the bivector part");
    check_scalar(ga::left_contraction(t, u), (t * u).scalar.value, "T _| U is the scalar part");
    check_vector(ga::left_contraction(t, q), (t * q).vector.e0, (t * q).vector.e1, (t * q).vector.e2, (t * q).vector.e3, "T _| I is the vector part");
    check_scalar(ga::left_contraction(q, q), 0.0, "I _| I is zero");

    // operator| is the left contraction.
    check_vector(v | b, ga::left_contraction(v, b).e0, ga::left_contraction(v, b).e1, ga::left_contraction(v, b).e2, ga::left_contraction(v, b).e3, "v | B");
    check_bivector(v | t, (v * t).bivector.e01, (v * t).bivector.e02, (v * t).bivector.e03, (v * t).bivector.e12, (v * t).bivector.e13, (v * t).bivector.e23, "v | T");
    check_trivector(v | q, (v * q).trivector.e012, (v * q).trivector.e013, (v * q).trivector.e023, (v * q).trivector.e123, "v | I");
    check_scalar(b | c, (b * c).scalar.value, "B | C");
    check_vector(b | t, (b * t).vector.e0, (b * t).vector.e1, (b * t).vector.e2, (b * t).vector.e3, "B | T");
    check_bivector(b | q, (b * q).bivector.e01, (b * q).bivector.e02, (b * q).bivector.e03, (b * q).bivector.e12, (b * q).bivector.e13, (b * q).bivector.e23, "B | I");
    check_scalar(t | u, (t * u).scalar.value, "T | U");
    check_vector(t | q, (t * q).vector.e0, (t * q).vector.e1, (t * q).vector.e2, (t * q).vector.e3, "T | I");
    check_scalar(q | q, 0.0, "I | I");

    // Right contraction mirrors it.
    check_scalar(ga::right_contraction(v, w), (v * w).scalar.value, "v |_ w");
    check_vector(ga::right_contraction(b, v), (b * v).vector.e0, (b * v).vector.e1, (b * v).vector.e2, (b * v).vector.e3, "B |_ v is the vector part");
    check_bivector(ga::right_contraction(t, v), (t * v).bivector.e01, (t * v).bivector.e02, (t * v).bivector.e03,
                   (t * v).bivector.e12, (t * v).bivector.e13, (t * v).bivector.e23, "T |_ v is the bivector part");
    check_trivector(ga::right_contraction(q, v), (q * v).trivector.e012, (q * v).trivector.e013, (q * v).trivector.e023,
                    (q * v).trivector.e123, "I |_ v is the trivector part");
    check_scalar(ga::right_contraction(b, c), (b * c).scalar.value, "B |_ C");
    check_vector(ga::right_contraction(t, b), (t * b).vector.e0, (t * b).vector.e1, (t * b).vector.e2, (t * b).vector.e3, "T |_ B");
    check_bivector(ga::right_contraction(q, b), (q * b).bivector.e01, (q * b).bivector.e02, (q * b).bivector.e03,
                   (q * b).bivector.e12, (q * b).bivector.e13, (q * b).bivector.e23, "I |_ B");
    check_scalar(ga::right_contraction(t, u), (t * u).scalar.value, "T |_ U");
    check_vector(ga::right_contraction(q, t), (q * t).vector.e0, (q * t).vector.e1, (q * t).vector.e2, (q * t).vector.e3, "I |_ T");
    check_scalar(ga::right_contraction(q, q), 0.0, "I |_ I");

    // The relations between the two: a _| B = -B |_ a for a vector and a bivector,
    // a _| T = T |_ a for a vector and a trivector.
    check_vector(ga::left_contraction(v, b), -ga::right_contraction(b, v).e0, -ga::right_contraction(b, v).e1,
                 -ga::right_contraction(b, v).e2, -ga::right_contraction(b, v).e3, "v _| B == -(B |_ v)");
    check_bivector(ga::left_contraction(v, t), ga::right_contraction(t, v).e01, ga::right_contraction(t, v).e02, ga::right_contraction(t, v).e03,
                   ga::right_contraction(t, v).e12, ga::right_contraction(t, v).e13, ga::right_contraction(t, v).e23, "v _| T == T |_ v");

    // scalar_product is the grade 0 part of anything.
    check_scalar(ga::scalar_product(v, w), (v | w).value, "scalar_product(v, w)");
    check_scalar(ga::scalar_product(b, c), (b * c).scalar.value, "scalar_product(B, C)");
    check_scalar(ga::scalar_product(v, b), 0.0, "scalar_product(v, B) is zero");

    // Where the product is already a pure Scalar, scalar_product must still
    // work. These did not compile until grade0 gained a Scalar overload.
    check_scalar(ga::scalar_product(Quadvector<double>(2.0), Quadvector<double>(3.0)), 0.0,
                 "scalar_product of two quadvectors, e0123^2 = 0");
    check_scalar(ga::scalar_product(Scalar<double>(2.0), Scalar<double>(5.0)), 10.0,
                 "scalar_product of two scalars");
    check_scalar(ga::grade0(Scalar<double>(7.0)), 7.0, "grade0(Scalar) is the identity");

    // Geometry: the plane x = 0 contracted into the point (1, 2, 3) is the line
    // through the point along x -- the same direction as e23, offset.
    const Bivector<double> perp = test_plane(1, 0, 0, 0) | test_point(1, 2, 3);
    check_close(perp.e12, 0.0, "x=0 _| P: no e12");
    check_close(perp.e13, 0.0, "x=0 _| P: no e13");
    check(perp.e23 != 0.0, "x=0 _| P: runs along x");
    check((ga::join(test_point(1, 2, 3), perp)).e1 == 0.0 && ga::norm(ga::join(test_point(1, 2, 3), perp)).value == 0.0,
          "x=0 _| P passes through P");
}

} // namespace

void run_grades_and_contractions()
{
    test_grade_projection();
    test_contractions();
}

} // namespace cctest
