// Cl(2,0) tests -- complex_isomorphism: every rotor operation checked against
// std::complex, which is an independent implementation of Cl(2,0)'s even
// subalgebra.
//
// One category of the Cl2 suite; main.cpp runs them all, and support.hpp
// holds what more than one file needs.
//
// This file -- and only this file -- includes <complex>. The even subalgebra of
// Cl(2,0) is the complex numbers, so std::complex is an independent reference
// implementation of every rotor operation, and test_complex_isomorphism checks
// the library against it. No header may take that dependency.

#include "support.hpp"

#include <complex>

namespace cctest {
namespace {

// ---------------------------------------------------------------------------
// The even subalgebra is the complex numbers
// ---------------------------------------------------------------------------

void test_complex_isomorphism()
{
    section("complex isomorphism");

    // std::complex is an independent implementation of the even subalgebra, so
    // every rotor operation is checked against it. The alias makes the reading
    // explicit.
    using C = std::complex<double>;
    const auto as_complex = [](const Rotor<double>& r) { return C(r.scalar.value, r.bivector.xy); };

    const ga::Complex<double> z1(Scalar<double>(1.5), Bivector<double>(-0.5));
    const ga::Complex<double> z2(Scalar<double>(-0.25), Bivector<double>(2.0));
    const C c1 = as_complex(z1), c2 = as_complex(z2);

    const C product = c1 * c2;
    check_rotor(z1 * z2, product.real(), product.imag(), "rotor product == complex product");

    const C quotient = C(1.0, 0.0) / c1;
    check_rotor(ga::inverse(z1), quotient.real(), quotient.imag(), "inverse == 1/z");

    check_close(ga::norm(z1).value, std::abs(c1), "norm == |z|");
    check_rotor(ga::reverse(z1), std::conj(c1).real(), std::conj(c1).imag(), "reverse == conj(z)");
    check_close(ga::log(z1).xy, std::arg(c1), "log == arg(z)");

    const C e = std::exp(C(0.0, 0.7));
    check_rotor(ga::exp(Bivector<double>(0.7)), e.real(), e.imag(), "exp(theta e12) == e^(i theta)");

    // And a rotation in the plane is "multiply by e^(i theta)" -- here as R^2.
    const C point(2.0, -1.5);
    const C turned = std::exp(C(0.0, 1.1)) * point;
    const Vector<double> rotated = ga::rotate(Vector<double>(2.0, -1.5), ga::rotor_from_angle(1.1));
    check_close(rotated.x, turned.real(), "rotate == multiply by e^(i theta) (x)");
    check_close(rotated.y, turned.imag(), "rotate == multiply by e^(i theta) (y)");
}

} // namespace

void run_complex_isomorphism()
{
    test_complex_isomorphism();
}

} // namespace cctest
