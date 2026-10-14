// Cl(3,0,1) tests -- arithmetic: compound assignment on all nine types, and
// mixed-grade sums and differences, including the explicit narrowings
// to_motor, to_rotor and to_translator.
//
// One category of the PGA suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

void test_compound_assignment()
{
    section("compound assignment");

    { Scalar<double> s(2.0); s += Scalar<double>(1.0); s -= Scalar<double>(0.5); s *= Scalar<double>(4.0); s /= Scalar<double>(2.0);
      check_scalar(s, 5.0, "Scalar += -= *= /="); }
    { Vector<double> v(1, 2, 3, 4); v += Vector<double>(1, 1, 1, 1); v -= Vector<double>(0, 1, 0, 1); v *= 2.0; v /= Scalar<double>(2.0);
      check_vector(v, 2, 2, 4, 4, "Vector += -= *= /="); }
    { Bivector<double> b(1, 2, 3, 4, 5, 6); b += b; b -= Bivector<double>(1, 1, 1, 1, 1, 1); b *= Scalar<double>(2.0); b /= 2.0;
      check_bivector(b, 1, 3, 5, 7, 9, 11, "Bivector += -= *= /="); }
    { Trivector<double> t(1, 2, 3, 4); t += t; t -= Trivector<double>(1, 1, 1, 1); t *= 2.0; t /= Scalar<double>(2.0);
      check_trivector(t, 1, 3, 5, 7, "Trivector += -= *= /="); }
    { Quadvector<double> q(3.0); q += q; q -= Quadvector<double>(1.0); q *= 2.0; q /= Scalar<double>(2.0);
      check_quadvector(q, 5.0, "Quadvector += -= *= /="); }
    { Multivector<double> m = make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);
      m += m; m -= make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16); m *= Scalar<double>(2.0); m /= Scalar<double>(2.0);
      check_mv_close(m, make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), "Multivector += -= *= /="); }
    { Rotor<double> r(Scalar<double>(1.0), 2.0, 3.0, 4.0); r += r; r -= Rotor<double>(Scalar<double>(1.0), 1, 1, 1); r *= Scalar<double>(2.0); r /= Scalar<double>(2.0);
      check_rotor(r, 1, 3, 5, 7, "Rotor += -= *= /="); }
    { Translator<double> t(Scalar<double>(1.0), 2.0, 3.0, 4.0); t += t; t -= Translator<double>(Scalar<double>(1.0), 1, 1, 1); t *= Scalar<double>(2.0); t /= Scalar<double>(2.0);
      check_translator(t, 1, 3, 5, 7, "Translator += -= *= /="); }
    { Motor<double> m(Scalar<double>(1.0), Bivector<double>(1, 1, 1, 1, 1, 1), Quadvector<double>(1.0));
      m += m; m -= Motor<double>(Scalar<double>(1.0), Bivector<double>(), Quadvector<double>()); m *= Scalar<double>(2.0); m /= Scalar<double>(2.0);
      check_scalar(m.scalar, 1.0, "Motor += -= *= /= scalar");
      check_bivector(m.bivector, 2, 2, 2, 2, 2, 2, "Motor += -= *= /= bivector");
      check_quadvector(m.quadvector, 2.0, "Motor += -= *= /= quadvector"); }
}

void test_mixed_grade_addition()
{
    section("mixed-grade addition");

    const Scalar<double> s(1.0);
    const Vector<double> v(2, 3, 4, 5);
    const Bivector<double> b(6, 7, 8, 9, 10, 11);
    const Trivector<double> t(12, 13, 14, 15);
    const Quadvector<double> q(16.0);
    const Multivector<double> all = make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);

    // A chain through every grade, in an arbitrary order, lands every component.
    check_mv_close(t + s + q + v + b, all, "s + v + b + t + q in any order");

    // Every distinct pair, both orders, against promote-and-add.
    check_mv_close(s + v, make_mv(1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), "Scalar + Vector");
    check_mv_close(v + s, s + v, "Vector + Scalar");
    check_mv_close(s + b, make_mv(1, 0, 0, 0, 0, 6, 7, 8, 9, 10, 11, 0, 0, 0, 0, 0), "Scalar + Bivector");
    check_mv_close(b + s, s + b, "Bivector + Scalar");
    check_mv_close(s + t, make_mv(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 13, 14, 15, 0), "Scalar + Trivector");
    check_mv_close(t + s, s + t, "Trivector + Scalar");
    check_mv_close(s + q, make_mv(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16), "Scalar + Quadvector");
    check_mv_close(q + s, s + q, "Quadvector + Scalar");
    check_mv_close(v + b, make_mv(0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 0, 0, 0, 0), "Vector + Bivector");
    check_mv_close(b + v, v + b, "Bivector + Vector");
    check_mv_close(v + t, make_mv(0, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 12, 13, 14, 15, 0), "Vector + Trivector");
    check_mv_close(t + v, v + t, "Trivector + Vector");
    check_mv_close(v + q, make_mv(0, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16), "Vector + Quadvector");
    check_mv_close(q + v, v + q, "Quadvector + Vector");
    check_mv_close(b + t, make_mv(0, 0, 0, 0, 0, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0), "Bivector + Trivector");
    check_mv_close(t + b, b + t, "Trivector + Bivector");
    check_mv_close(b + q, make_mv(0, 0, 0, 0, 0, 6, 7, 8, 9, 10, 11, 0, 0, 0, 0, 16), "Bivector + Quadvector");
    check_mv_close(q + b, b + q, "Quadvector + Bivector");
    check_mv_close(t + q, make_mv(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 13, 14, 15, 16), "Trivector + Quadvector");
    check_mv_close(q + t, t + q, "Quadvector + Trivector");

    // Multivector with each grade, both orders.
    check_mv_close(all + s, make_mv(2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), "Multivector + Scalar");
    check_mv_close(s + all, all + s, "Scalar + Multivector");
    check_mv_close(all + v, make_mv(1, 4, 6, 8, 10, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), "Multivector + Vector");
    check_mv_close(v + all, all + v, "Vector + Multivector");
    check_mv_close(all + b, make_mv(1, 2, 3, 4, 5, 12, 14, 16, 18, 20, 22, 12, 13, 14, 15, 16), "Multivector + Bivector");
    check_mv_close(b + all, all + b, "Bivector + Multivector");
    check_mv_close(all + t, make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 24, 26, 28, 30, 16), "Multivector + Trivector");
    check_mv_close(t + all, all + t, "Trivector + Multivector");
    check_mv_close(all + q, make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 32), "Multivector + Quadvector");
    check_mv_close(q + all, all + q, "Quadvector + Multivector");

    // Every mixed sum is a Multivector, even the ones that would fit a motor.
    static_assert(std::is_same<decltype(s + b), Multivector<double>>::value, "Scalar + Bivector widens to Multivector");
    static_assert(std::is_same<decltype(s + q), Multivector<double>>::value, "Scalar + Quadvector widens to Multivector");
    static_assert(std::is_same<decltype(v + t), Multivector<double>>::value, "Vector + Trivector widens to Multivector");
    static_assert(std::is_same<decltype(all + q), Multivector<double>>::value, "Multivector + Quadvector stays Multivector");

    // The even types widen into Multivector on the right.
    const Rotor<double> r(Scalar<double>(1.0), 2, 3, 4);
    const Translator<double> tr(Scalar<double>(1.0), 2, 3, 4);
    const Motor<double> mo(Scalar<double>(1.0), Bivector<double>(1, 1, 1, 1, 1, 1), Quadvector<double>(1.0));
    check_mv_close(all + r, make_mv(2, 2, 3, 4, 5, 6, 7, 8, 11, 13, 15, 12, 13, 14, 15, 16), "Multivector + Rotor");
    check_mv_close(all + tr, make_mv(2, 2, 3, 4, 5, 8, 10, 12, 9, 10, 11, 12, 13, 14, 15, 16), "Multivector + Translator");
    check_mv_close(all + mo, make_mv(2, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 12, 13, 14, 15, 17), "Multivector + Motor");

    // The named conversions.
    const Motor<double> rm = ga::to_motor(r);
    check_bivector(rm.bivector, 0, 0, 0, 2, 3, 4, "to_motor(Rotor)");
    const Motor<double> tm = ga::to_motor(tr);
    check_bivector(tm.bivector, 2, 3, 4, 0, 0, 0, "to_motor(Translator)");
    const Motor<double> am = ga::to_motor(all);
    check_scalar(am.scalar, 1.0, "to_motor(Multivector) scalar");
    check_bivector(am.bivector, 6, 7, 8, 9, 10, 11, "to_motor(Multivector) bivector");
    check_quadvector(am.quadvector, 16.0, "to_motor(Multivector) quadvector");
    check_rotor(ga::to_rotor(all), 1, 9, 10, 11, "to_rotor(Multivector)");
    check_rotor(ga::to_rotor(mo), 1, 1, 1, 1, "to_rotor(Motor)");
    check_translator(ga::to_translator(all), 1, 6, 7, 8, "to_translator(Multivector)");
    check_translator(ga::to_translator(mo), 1, 1, 1, 1, "to_translator(Motor)");
    check_mv_close(ga::to_multivector(r), make_mv(1, 0, 0, 0, 0, 0, 0, 0, 2, 3, 4, 0, 0, 0, 0, 0), "to_multivector(Rotor)");
    check_mv_close(ga::to_multivector(tr), make_mv(1, 0, 0, 0, 0, 2, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0), "to_multivector(Translator)");
    check_mv_close(ga::to_multivector(mo), make_mv(1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1), "to_multivector(Motor)");
}

void test_mixed_grade_subtraction()
{
    section("mixed-grade subtraction");

    const Scalar<double> s(1.0);
    const Vector<double> v(2, 3, 4, 5);
    const Bivector<double> b(6, 7, 8, 9, 10, 11);
    const Trivector<double> t(12, 13, 14, 15);
    const Quadvector<double> q(16.0);
    const Multivector<double> all = make_mv(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16);

    // a - b == a + (-b) for every pair, both orders.
    check_mv_close(s - v, s + (-v), "Scalar - Vector");     check_mv_close(v - s, v + (-s), "Vector - Scalar");
    check_mv_close(s - b, s + (-b), "Scalar - Bivector");   check_mv_close(b - s, b + (-s), "Bivector - Scalar");
    check_mv_close(s - t, s + (-t), "Scalar - Trivector");  check_mv_close(t - s, t + (-s), "Trivector - Scalar");
    check_mv_close(s - q, s + (-q), "Scalar - Quadvector"); check_mv_close(q - s, q + (-s), "Quadvector - Scalar");
    check_mv_close(v - b, v + (-b), "Vector - Bivector");   check_mv_close(b - v, b + (-v), "Bivector - Vector");
    check_mv_close(v - t, v + (-t), "Vector - Trivector");  check_mv_close(t - v, t + (-v), "Trivector - Vector");
    check_mv_close(v - q, v + (-q), "Vector - Quadvector"); check_mv_close(q - v, q + (-v), "Quadvector - Vector");
    check_mv_close(b - t, b + (-t), "Bivector - Trivector"); check_mv_close(t - b, t + (-b), "Trivector - Bivector");
    check_mv_close(b - q, b + (-q), "Bivector - Quadvector"); check_mv_close(q - b, q + (-b), "Quadvector - Bivector");
    check_mv_close(t - q, t + (-q), "Trivector - Quadvector"); check_mv_close(q - t, q + (-t), "Quadvector - Trivector");
    check_mv_close(all - s, all + (-s), "Multivector - Scalar");     check_mv_close(s - all, s + (-all), "Scalar - Multivector");
    check_mv_close(all - v, all + (-v), "Multivector - Vector");     check_mv_close(v - all, v + (-all), "Vector - Multivector");
    check_mv_close(all - b, all + (-b), "Multivector - Bivector");   check_mv_close(b - all, b + (-all), "Bivector - Multivector");
    check_mv_close(all - t, all + (-t), "Multivector - Trivector");  check_mv_close(t - all, t + (-all), "Trivector - Multivector");
    check_mv_close(all - q, all + (-q), "Multivector - Quadvector"); check_mv_close(q - all, q + (-all), "Quadvector - Multivector");

    // Subtracting a grade removes it.
    check_mv_close(all - v - b - t - q - s, Multivector<double>(), "removing every grade leaves zero");
    static_assert(std::is_same<decltype(v - q), Multivector<double>>::value, "Vector - Quadvector widens to Multivector");
}

} // namespace

void run_arithmetic()
{
    test_compound_assignment();
    test_mixed_grade_addition();
    test_mixed_grade_subtraction();
}

} // namespace cctest
