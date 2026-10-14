// Cl(3,0,1) tests -- comparison_and_printing: exact == and !=, approx_equal,
// and operator<<.
//
// One category of the PGA suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.

#include "support.hpp"

namespace cctest {
namespace {

void test_comparison()
{
    section("comparison");

    const Scalar<double> s(1.5); const Vector<double> v(1, 2, 3, 4); const Bivector<double> b(1, 2, 3, 4, 5, 6);
    const Trivector<double> t(1, 2, 3, 4); const Quadvector<double> q(5);
    const Multivector<double> m(s, v, b, t, q);
    const Rotor<double> r(s, 1, 2, 3); const Translator<double> tr(s, 1, 2, 3); const Motor<double> mo(s, b, q);

    check(s == Scalar<double>(1.5) && !(s != Scalar<double>(1.5)), "Scalar == and !=");
    check(v == Vector<double>(1, 2, 3, 4) && v != Vector<double>(1, 2, 3, 5), "Vector == and !=");
    check(b == Bivector<double>(1, 2, 3, 4, 5, 6) && b != Bivector<double>(1, 2, 3, 4, 5, 7), "Bivector == and !=");
    check(t == Trivector<double>(1, 2, 3, 4) && t != Trivector<double>(0, 2, 3, 4), "Trivector == and !=");
    check(q == Quadvector<double>(5) && q != Quadvector<double>(-5), "Quadvector == and !=");
    check(m == Multivector<double>(s, v, b, t, q) && m != Multivector<double>(s, v, b, t, Quadvector<double>()), "Multivector == and !=");
    check(r == Rotor<double>(s, 1, 2, 3) && r != Rotor<double>(s, 1, 2, 4), "Rotor == and !=");
    check(tr == Translator<double>(s, 1, 2, 3) && tr != Translator<double>(s, 1, 2, 4), "Translator == and !=");
    check(mo == Motor<double>(s, b, q) && mo != Motor<double>(s, b, Quadvector<double>()), "Motor == and !=");
    check(!(mo == -mo), "a motor and its negation compare unequal, as in Cl3");

    // approx_equal: default tolerance absorbs rounding, explicit tolerance is respected.
    const Vector<double> nudged(1, 2, 3, 4 + 1e-15);
    check(v != nudged, "exact == sees a 1e-15 difference");
    check(ga::approx_equal(v, nudged), "approx_equal absorbs it at the default tolerance");
    check(!ga::approx_equal(v, Vector<double>(1, 2, 3, 4.01)), "approx_equal rejects 0.01 at the default tolerance");
    check(ga::approx_equal(v, Vector<double>(1, 2, 3, 4.01), 0.1), "... but accepts it at 0.1");
    check(ga::approx_equal(s, Scalar<double>(1.5 + 1e-15)), "approx_equal(Scalar)");
    check(ga::approx_equal(b, b) && ga::approx_equal(t, t) && ga::approx_equal(q, q), "approx_equal on b, t, q");
    check(ga::approx_equal(m, m) && ga::approx_equal(r, r) && ga::approx_equal(tr, tr) && ga::approx_equal(mo, mo), "approx_equal on m, r, tr, mo");
    check(!ga::approx_equal(mo, -mo), "approx_equal(motor, -motor) is false");

    // Streams delegate to to_string.
    std::ostringstream os;
    os << s << " " << v << " " << b << " " << t << " " << q << " " << m << " " << r << " " << tr << " " << mo;
    check(os.str().find("e0123") != std::string::npos, "operator<< prints the quadvector label");
    check(os.str() == s.to_string() + " " + v.to_string() + " " + b.to_string() + " " + t.to_string() + " " + q.to_string()
          + " " + m.to_string() + " " + r.to_string() + " " + tr.to_string() + " " + mo.to_string(),
          "operator<< is to_string for every type");
}

} // namespace

void run_comparison_and_printing()
{
    test_comparison();
}

} // namespace cctest
