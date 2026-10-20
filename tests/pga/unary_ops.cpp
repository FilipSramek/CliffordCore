// Cl(3,0,1) tests -- unary_ops: norm and ideal_norm, normalize, the
// involutions, inverse, and dual/undual, for every type.
//
// One category of the PGA suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

void test_norms()
{
    section("norm");

    // Euclidean norms ignore every e0 component; ideal norms see only those.
    const Vector<double> v(7, 3, 4, 0);   // e0 = 7, normal (3, 4, 0)
    check_scalar(ga::norm(v), 5.0, "norm(plane) is the length of the normal");
    check_scalar(ga::squared_norm(v), 25.0, "squared_norm(plane)");
    check_scalar(ga::ideal_norm(v), 7.0, "ideal_norm(plane) is |e0|");
    check_scalar(ga::norm(Vector<double>(2, 0, 0, 0)), 0.0, "the ideal plane has zero norm");

    const Bivector<double> b(1, 2, 2, 0, 3, 4);
    check_scalar(ga::norm(b), 5.0, "norm(line) is the length of the direction");
    check_scalar(ga::ideal_norm(b), 3.0, "ideal_norm(line) is the size of the moment");
    check_scalar(ga::norm(Bivector<double>(1, 2, 3, 0, 0, 0)), 0.0, "an ideal line has zero norm");

    const Trivector<double> t(0, 3, 4, -2);
    check_scalar(ga::norm(t), 2.0, "norm(point) is the absolute weight");
    check_scalar(ga::squared_norm(t), 4.0, "squared_norm(point)");
    check_scalar(ga::ideal_norm(t), 5.0, "ideal_norm(point) is distance times weight");
    check_scalar(ga::norm(test_point(3, 4, 0)), 1.0, "a unit point has norm 1");
    check_scalar(ga::ideal_norm(test_point(3, 4, 0)), 5.0, "... and ideal norm equal to its distance from the origin");

    const Quadvector<double> q(-3.0);
    check_scalar(ga::norm(q), 0.0, "norm(quadvector) is always zero");
    check_scalar(ga::squared_norm(q), 0.0, "squared_norm(quadvector) is always zero");
    check_scalar(ga::ideal_norm(q), 3.0, "ideal_norm(quadvector) is |e0123|");

    const Scalar<double> s(-2.5);
    check_scalar(ga::norm(s), 2.5, "norm(scalar)");
    check_scalar(ga::ideal_norm(s), 0.0, "ideal_norm(scalar) is zero");

    // Multivector: the non-e0 squares, and the scalar part of m * reverse(m).
    const Multivector<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
    const double euclid = 1 + 9 + 16 + 25 + 81 + 100 + 121 + 225;
    const double ideal = 4 + 36 + 49 + 64 + 144 + 169 + 196 + 256;
    check_scalar(ga::squared_norm(m), euclid, "squared_norm(multivector) sums the non-e0 squares");
    check_scalar(ga::norm(m), std::sqrt(euclid), "norm(multivector)");
    check_scalar(ga::ideal_norm(m), std::sqrt(ideal), "ideal_norm(multivector) sums the e0 squares");
    check_close((m * ga::reverse(m)).scalar.value, euclid, "squared_norm is the scalar part of m * reverse(m)");

    const Rotor<double> r(Scalar<double>(2.0), 3, 6, 0);
    check_scalar(ga::norm(r), 7.0, "norm(rotor)");
    check_scalar(ga::ideal_norm(r), 0.0, "ideal_norm(rotor) is zero");
    const Translator<double> tr(Scalar<double>(-3.0), 1, 2, 2);
    check_scalar(ga::norm(tr), 3.0, "norm(translator) is |scalar|");
    check_scalar(ga::ideal_norm(tr), 3.0, "ideal_norm(translator)");
    const Motor<double> mo(Scalar<double>(2.0), Bivector<double>(1, 2, 2, 3, 6, 0), Quadvector<double>(4.0));
    check_scalar(ga::norm(mo), 7.0, "norm(motor) ignores the ideal part and e0123");
    check_scalar(ga::ideal_norm(mo), 5.0, "ideal_norm(motor)");
    check_close((mo * ga::reverse(mo)).scalar.value, 49.0, "squared_norm(motor) is the scalar part of the study number");
}

void test_normalize()
{
    section("normalize");

    check_scalar(ga::normalize(Scalar<double>(-4.0)), -1.0, "normalize(scalar) is its sign");
    check_scalar(ga::norm(ga::normalize(Vector<double>(7, 3, 4, 0))), 1.0, "normalize(plane) has a unit normal");
    check_vector(ga::normalize(Vector<double>(10, 0, 0, 5)), 2, 0, 0, 1, "normalize(plane) scales the offset along");
    check_scalar(ga::norm(ga::normalize(Bivector<double>(1, 2, 2, 0, 3, 4))), 1.0, "normalize(line) has a unit direction");
    check_trivector(ga::normalize(Trivector<double>(-6, 4, -2, 2)), -3, 2, -1, 1, "normalize(point) has weight 1");
    check_scalar(ga::norm(ga::normalize(make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16))), 1.0,
                 "normalize(multivector) has unit norm");
    check_scalar(ga::norm(ga::normalize(Rotor<double>(Scalar<double>(2.0), 3, 6, 0))), 1.0, "normalize(rotor) has unit norm");
    check_translator(ga::normalize(Translator<double>(Scalar<double>(2.0), 1, -2, 4)), 1, 0.5, -1, 2,
                     "normalize(translator) has scalar 1");

    // A drifted motor: dividing by the norm alone would leave m * reverse(m)
    // with an e0123 residue. The study-number normalisation removes it.
    const Motor<double> drifted(Scalar<double>(2.0), Bivector<double>(1, 2, 2, 3, 6, 0), Quadvector<double>(4.0));
    const Motor<double> unit = ga::normalize(drifted);
    const Motor<double> check_ = unit * ga::reverse(unit);
    check_scalar(check_.scalar, 1.0, "normalize(motor): scalar part of m * reverse(m) is 1");
    check_quadvector(check_.quadvector, 0.0, "normalize(motor): e0123 part of m * reverse(m) is 0");
    check_bivector(check_.bivector, 0, 0, 0, 0, 0, 0, "normalize(motor): m * reverse(m) has no bivector part");
    // The naive scaling really would have been wrong for this motor.
    const Motor<double> naive = drifted / ga::norm(drifted);
    check(std::fabs((naive * ga::reverse(naive)).quadvector.e0123) > 1e-3, "dividing by the norm alone leaves an e0123 residue");
    // And an already-unit motor is left alone.
    check_mv_close(ga::to_multivector(ga::normalize(unit)), ga::to_multivector(unit), "normalize is idempotent on a unit motor");
    // A pure rotor and a pure translator normalise the same way through Motor.
    const Motor<double> fromRotor = ga::normalize(Motor<double>(Rotor<double>(Scalar<double>(2.0), 3, 6, 0)));
    check_close(ga::norm(fromRotor).value, 1.0, "normalize(Motor(rotor)) has unit norm");
    check_quadvector(fromRotor.quadvector, 0.0, "normalize(Motor(rotor)) gains no e0123");
}

#ifdef NDEBUG
// Release builds only: with the assert compiled out, normalize() hands an
// element with no Euclidean size -- zero, or ideal -- back unchanged. See the
// same test in tests/cl3/unary_ops.cpp, and audit items A1 and A4 in
// docs/safety-audit.md.
void test_normalize_zero_release()
{
    section("normalize of zero or ideal, release build");

    check(ga::normalize(Scalar<double>()) == Scalar<double>(), "normalize(Scalar 0) returns it unchanged");
    const Vector<double> idealPlane(5.0, 0.0, 0.0, 0.0);
    check(ga::normalize(idealPlane) == idealPlane, "normalize(ideal plane) returns it unchanged");
    const Bivector<double> idealLine(1.0, 2.0, 3.0, 0.0, 0.0, 0.0);
    check(ga::normalize(idealLine) == idealLine, "normalize(ideal line) returns it unchanged");
    const Trivector<double> idealPoint(1.0, 2.0, 3.0, 0.0);
    check(ga::normalize(idealPoint) == idealPoint, "normalize(ideal point) returns it unchanged");
    check(ga::normalize(Multivector<double>()) == Multivector<double>(), "normalize(Multivector 0) returns it unchanged");
    const Rotor<double> zeroRotor(Scalar<double>(0.0), 0.0, 0.0, 0.0);
    check(ga::normalize(zeroRotor) == zeroRotor, "normalize(Rotor 0) returns it unchanged");
    const Translator<double> idealTranslator(Scalar<double>(0.0), 1.0, 2.0, 3.0);
    check(ga::normalize(idealTranslator) == idealTranslator, "normalize(ideal translator) returns it unchanged");
    const Motor<double> idealMotor(Scalar<double>(0.0), idealLine, Quadvector<double>(1.0));
    check(ga::normalize(idealMotor) == idealMotor, "normalize(ideal motor) returns it unchanged");
    // A plane that is NOT ideal takes the same branch once its squared norm
    // underflows -- audit item A4. (A point would not: its norm is |e123|,
    // with no square to underflow.)
    const Vector<double> tiny(0.0, 1e-200, 0.0, 0.0);
    check(ga::normalize(tiny) == tiny, "normalize(plane with normal 1e-200) returns it unchanged, not unit");
}
#endif

void test_involutions()
{
    section("involutions");

    // Sign per grade: reverse + + - - +, involute + - + - +, conjugate + - - + +.
    const Multivector<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
    check_mv_close(ga::reverse(m), make_mv(1, 2, 3, 4, 5, -6, -7, -8, -9, -10, -11, -12, -13, -14, -15, 16),
                   "reverse negates grades 2 and 3");
    check_mv_close(ga::involute(m), make_mv(1, -2, -3, -4, -5, 6, 7, 8, 9, 10, 11, -12, -13, -14, -15, 16),
                   "involute negates grades 1 and 3");
    check_mv_close(ga::conjugate(m), make_mv(1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11, 12, 13, 14, 15, 16),
                   "conjugate negates grades 1 and 2");
    check_mv_close(ga::conjugate(m), ga::reverse(ga::involute(m)), "conjugate == reverse of involute");
    check_mv_close(ga::conjugate(m), ga::involute(ga::reverse(m)), "conjugate == involute of reverse");
    check_mv_close(ga::reverse(ga::reverse(m)), m, "reverse is an involution");

    // Anti-automorphism: reverse(ab) == reverse(b) reverse(a).
    const Multivector<double> n = make_mv(-0.5, 1, 2, -3, 4, -1, 0.5, 2, 3, -2, 1, 0.5, 1.5, -1, 1, -2);
    check_mv_close(ga::reverse(m * n), ga::reverse(n) * ga::reverse(m), "reverse(ab) == reverse(b) reverse(a)", 1e-9);
    check_mv_close(ga::involute(m * n), ga::involute(m) * ga::involute(n), "involute(ab) == involute(a) involute(b)", 1e-9);
    check_mv_close(ga::conjugate(m * n), ga::conjugate(n) * ga::conjugate(m), "conjugate(ab) == conjugate(b) conjugate(a)", 1e-9);

    // Per-type overloads agree with the multivector ones.
    const Scalar<double> s(3); const Vector<double> v(1, 2, 3, 4); const Bivector<double> b(1, 2, 3, 4, 5, 6);
    const Trivector<double> t(1, 2, 3, 4); const Quadvector<double> q(5);
    check_scalar(ga::reverse(s), 3, "reverse(scalar)");      check_scalar(ga::involute(s), 3, "involute(scalar)");    check_scalar(ga::conjugate(s), 3, "conjugate(scalar)");
    check_vector(ga::reverse(v), 1, 2, 3, 4, "reverse(vector)"); check_vector(ga::involute(v), -1, -2, -3, -4, "involute(vector)"); check_vector(ga::conjugate(v), -1, -2, -3, -4, "conjugate(vector)");
    check_bivector(ga::reverse(b), -1, -2, -3, -4, -5, -6, "reverse(bivector)"); check_bivector(ga::involute(b), 1, 2, 3, 4, 5, 6, "involute(bivector)"); check_bivector(ga::conjugate(b), -1, -2, -3, -4, -5, -6, "conjugate(bivector)");
    check_trivector(ga::reverse(t), -1, -2, -3, -4, "reverse(trivector)"); check_trivector(ga::involute(t), -1, -2, -3, -4, "involute(trivector)"); check_trivector(ga::conjugate(t), 1, 2, 3, 4, "conjugate(trivector)");
    check_quadvector(ga::reverse(q), 5, "reverse(quadvector)"); check_quadvector(ga::involute(q), 5, "involute(quadvector)"); check_quadvector(ga::conjugate(q), 5, "conjugate(quadvector)");

    const Rotor<double> r(Scalar<double>(1), 2, 3, 4);
    check_rotor(ga::reverse(r), 1, -2, -3, -4, "reverse(rotor)");
    check_rotor(ga::involute(r), 1, 2, 3, 4, "involute(rotor) is the identity");
    check_rotor(ga::conjugate(r), 1, -2, -3, -4, "conjugate(rotor) == reverse(rotor)");
    const Translator<double> tr(Scalar<double>(1), 2, 3, 4);
    check_translator(ga::reverse(tr), 1, -2, -3, -4, "reverse(translator)");
    check_translator(ga::involute(tr), 1, 2, 3, 4, "involute(translator) is the identity");
    check_translator(ga::conjugate(tr), 1, -2, -3, -4, "conjugate(translator) == reverse(translator)");
    const Motor<double> mo(Scalar<double>(1), Bivector<double>(2, 3, 4, 5, 6, 7), Quadvector<double>(8));
    check_mv_close(ga::to_multivector(ga::reverse(mo)), ga::reverse(ga::to_multivector(mo)), "reverse(motor) matches reverse(multivector)");
    check_mv_close(ga::to_multivector(ga::involute(mo)), ga::involute(ga::to_multivector(mo)), "involute(motor) matches involute(multivector)");
    check_mv_close(ga::to_multivector(ga::conjugate(mo)), ga::conjugate(ga::to_multivector(mo)), "conjugate(motor) matches conjugate(multivector)");
}

void test_inverse()
{
    section("inverse");

    const Multivector<double> one = basis_blade(0);

    check_scalar(ga::inverse(Scalar<double>(4.0)), 0.25, "inverse(scalar)");

    const Vector<double> v(7, 3, 4, 0);
    check_mv_close(v * ga::inverse(v), one, "plane * inverse(plane) == 1");
    check_mv_close(ga::inverse(v) * v, one, "inverse(plane) * plane == 1");
    check_vector(ga::inverse(Vector<double>(5, 0, 0, 0)), 0, 0, 0, 0, "the ideal plane has no inverse: zero is returned");

    // A general bivector is not a line and carries an e0123 term in B * reverse(B).
    const Bivector<double> b(1, 2, 2, 0, 3, 4);
    check_mv_close(b * ga::inverse(b), one, "bivector * inverse(bivector) == 1");
    check_mv_close(ga::inverse(b) * b, one, "inverse(bivector) * bivector == 1");
    const Bivector<double> zAxis(0, 0, 0, 1, 0, 0);
    check_bivector(ga::inverse(zAxis), 0, 0, 0, -1, 0, 0, "a unit line inverts to its negation");
    check_bivector(ga::inverse(Bivector<double>(1, 2, 3, 0, 0, 0)), 0, 0, 0, 0, 0, 0, "an ideal line has no inverse: zero is returned");

    const Trivector<double> t(-6, 4, -2, 2);
    check_mv_close(t * ga::inverse(t), one, "point * inverse(point) == 1");
    check_trivector(ga::inverse(test_point(1, 2, 3)), 3, -2, 1, -1, "a unit point inverts to its negation");
    check_trivector(ga::inverse(Trivector<double>(1, 2, 3, 0)), 0, 0, 0, 0, "an ideal point has no inverse: zero is returned");

    // The full 16-component case, through the Hitzer-Sangwine formula.
    const Multivector<double> m = make_mv(1, -2, 3, 0.5, -1, 2, -3, 1, 0.25, -0.5, 4, 1, -1, 2, -2, 0.75);
    check_mv_close(m * ga::inverse(m), one, "multivector * inverse(multivector) == 1", 1e-9);
    check_mv_close(ga::inverse(m) * m, one, "inverse(multivector) * multivector == 1", 1e-9);
    check_mv_close(ga::inverse(basis_blade(1)), Multivector<double>(), "e0 has no inverse: zero is returned");
    check_mv_close(ga::inverse(basis_blade(15)), Multivector<double>(), "e0123 has no inverse: zero is returned");
    // The formula reduces to the simple cases.
    check_mv_close(ga::inverse(ga::detail::promote(v)), ga::detail::promote(ga::inverse(v)), "inverse(Multivector) agrees with inverse(Vector)");
    check_mv_close(ga::inverse(ga::detail::promote(b)), ga::detail::promote(ga::inverse(b)), "inverse(Multivector) agrees with inverse(Bivector)");

    const Rotor<double> r(Scalar<double>(2.0), 3, 6, 0);
    check_mv_close(ga::to_multivector(r * ga::inverse(r)), one, "rotor * inverse(rotor) == 1");
    const Translator<double> tr(Scalar<double>(-3.0), 1, 2, 2);
    check_mv_close(ga::to_multivector(tr * ga::inverse(tr)), one, "translator * inverse(translator) == 1");
    check_translator(ga::inverse(Translator<double>(Scalar<double>(1.0), 1, 2, 3)), 1, -1, -2, -3,
                     "a normalised translator inverts to its reverse");

    // An UNNORMALISED motor with a non-zero study part: reverse/|m|^2 is wrong here.
    const Motor<double> mo(Scalar<double>(2.0), Bivector<double>(1, 2, 2, 3, 6, 0), Quadvector<double>(4.0));
    check_mv_close(ga::to_multivector(mo * ga::inverse(mo)), one, "motor * inverse(motor) == 1");
    check_mv_close(ga::to_multivector(ga::inverse(mo) * mo), one, "inverse(motor) * motor == 1");
    const Motor<double> naive = ga::reverse(mo) / ga::squared_norm(mo);
    check(mv_difference(ga::to_multivector(mo * naive), one) > 1e-3, "reverse(m)/|m|^2 is NOT the inverse of a drifted motor");
    check_mv_close(ga::to_multivector(ga::inverse(mo)), ga::inverse(ga::to_multivector(mo)), "inverse(Motor) agrees with inverse(Multivector)", 1e-9);
}

void test_dual()
{
    section("dual");

    // The full table from docs/pga.md: blade index -> (sign, blade index).
    const int expectedBlade[16] = {15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
    const double expectedSign[16] = {+1, +1, -1, +1, -1, +1, -1, +1, +1, -1, +1, +1, -1, +1, -1, +1};
    for (int i = 0; i < 16; ++i) {
        const Multivector<double> d = ga::dual(basis_blade(i));
        check_mv_close(d, expectedSign[i] * basis_blade(expectedBlade[i]),
                       std::string("dual(") + kBladeNames[i] + ") == " + (expectedSign[i] > 0 ? "+" : "-") + kBladeNames[expectedBlade[i]]);
        // a ^ dual(a) == +e0123 is the defining property. The wedge of two
        // multivectors is not an operator, so take the grade 4 part of the
        // product, which is the wedge for two blades of complementary grade.
        check_quadvector((basis_blade(i) * d).quadvector, 1.0, std::string(kBladeNames[i]) + " ^ dual == +e0123");
    }

    // Double dual: +x on even grades, -x on odd; undual undoes exactly.
    const Multivector<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
    check_mv_close(ga::dual(ga::dual(m)), make_mv(1, -2, -3, -4, -5, 6, 7, 8, 9, 10, 11, -12, -13, -14, -15, 16),
                   "dual(dual(x)) flips the odd grades");
    check_mv_close(ga::undual(ga::dual(m)), m, "undual(dual(x)) == x");
    check_mv_close(ga::dual(ga::undual(m)), m, "dual(undual(x)) == x");

    // Per-type overloads agree with the multivector one.
    const Scalar<double> s(3); const Vector<double> v(1, 2, 3, 4); const Bivector<double> b(1, 2, 3, 4, 5, 6);
    const Trivector<double> t(1, 2, 3, 4); const Quadvector<double> q(5);
    check_mv_close(ga::detail::promote(ga::dual(s)), ga::dual(ga::detail::promote(s)), "dual(scalar)");
    check_mv_close(ga::detail::promote(ga::dual(v)), ga::dual(ga::detail::promote(v)), "dual(vector)");
    check_mv_close(ga::detail::promote(ga::dual(b)), ga::dual(ga::detail::promote(b)), "dual(bivector)");
    check_mv_close(ga::detail::promote(ga::dual(t)), ga::dual(ga::detail::promote(t)), "dual(trivector)");
    check_mv_close(ga::detail::promote(ga::dual(q)), ga::dual(ga::detail::promote(q)), "dual(quadvector)");
    check_mv_close(ga::detail::promote(ga::undual(s)), ga::undual(ga::detail::promote(s)), "undual(scalar)");
    check_mv_close(ga::detail::promote(ga::undual(v)), ga::undual(ga::detail::promote(v)), "undual(vector)");
    check_mv_close(ga::detail::promote(ga::undual(b)), ga::undual(ga::detail::promote(b)), "undual(bivector)");
    check_mv_close(ga::detail::promote(ga::undual(t)), ga::undual(ga::detail::promote(t)), "undual(trivector)");
    check_mv_close(ga::detail::promote(ga::undual(q)), ga::undual(ga::detail::promote(q)), "undual(quadvector)");

    // Geometry: the dual of a plane through the origin is its normal direction.
    check_trivector(ga::dual(test_plane(1, 0, 0, 0)), 0, 0, -1, 0, "dual(x=0) is the ideal point -e023 (the x direction)");
    check_vector(ga::dual(test_point(0, 0, 0)), -1, 0, 0, 0, "dual(origin) is the ideal plane -e0");
    check_scalar(ga::ideal_norm(ga::dual(v)), ga::norm(v).value, "dual swaps the Euclidean and ideal norms");
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
