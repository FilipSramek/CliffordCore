#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include "../include/vector3.hpp"
#include "../include/operations/dot_product.hpp"
#include "../include/operations/wedge_product.hpp"
#include "../include/operations/geometric_product.hpp"
#include "../include/operations/norm.hpp"
#include "../include/operations/inverse.hpp"
#include "../include/operations/reverse.hpp"
#include "../include/operations/dual.hpp"

using CliffordCore::Bivector3;
using CliffordCore::Vector3;

namespace {

template <typename T>
void expect_close(T actual, T expected, T tolerance, const std::string& message)
{
    if (std::fabs(actual - expected) > tolerance) {
        std::cerr << "FAIL: " << message << " (expected " << expected << ", got " << actual << ")\n";
        std::exit(1);
    }
}

template <typename T>
void expect_true(bool condition, const std::string& message)
{
    if (!condition) {
        std::cerr << "FAIL: " << message << "\n";
        std::exit(1);
    }
}

} // namespace

int main()
{
    Vector3<double> a(1.0, 2.0, 3.0);
    Vector3<double> b(4.0, 5.0, 6.0);

    auto dot = CliffordCore::dot_product(a, b);
    expect_close(dot.value, 32.0, 1e-12, "dot product");

    auto wedge = CliffordCore::wedge_product(a, b);
    expect_close(wedge.xy, -3.0, 1e-12, "wedge xy");
    expect_close(wedge.xz, -6.0, 1e-12, "wedge xz");
    expect_close(wedge.yz, -3.0, 1e-12, "wedge yz");

    auto product = CliffordCore::geometric_product(a, b);
    expect_close(product.scalar.value, 32.0, 1e-12, "geometric scalar part");
    expect_close(product.bivector.xy, -3.0, 1e-12, "geometric bivector xy");
    expect_close(product.bivector.xz, -6.0, 1e-12, "geometric bivector xz");
    expect_close(product.bivector.yz, -3.0, 1e-12, "geometric bivector yz");

    auto norm_a = CliffordCore::norm(a);
    expect_close(norm_a.value, std::sqrt(14.0), 1e-12, "norm");

    auto inv_a = CliffordCore::inverse(a);
    expect_close(inv_a.x, 1.0 / 14.0, 1e-12, "inverse x");
    expect_close(inv_a.y, 2.0 / 14.0, 1e-12, "inverse y");
    expect_close(inv_a.z, 3.0 / 14.0, 1e-12, "inverse z");

    auto reversed_bivector = CliffordCore::reverse(Bivector3<double>(1.0, 2.0, 3.0));
    expect_close(reversed_bivector.xy, -1.0, 1e-12, "reverse bivector xy");
    expect_close(reversed_bivector.xz, -2.0, 1e-12, "reverse bivector xz");
    expect_close(reversed_bivector.yz, -3.0, 1e-12, "reverse bivector yz");

    auto dual_vec = CliffordCore::dual(a);
    expect_close(dual_vec.xy, 1.0, 1e-12, "dual vector xy");
    expect_close(dual_vec.xz, 2.0, 1e-12, "dual vector xz");
    expect_close(dual_vec.yz, 3.0, 1e-12, "dual vector yz");

    std::cout << "All CliffordCore tests passed.\n";
    return 0;
}
