#include <iostream>
#include "vector3.hpp"
#include "operations/dot_product.hpp"
#include "operations/wedge_product.hpp"
#include "operations/geometric_product.hpp"

int main() {
    CliffordCore::Vector3<double> a(1.0, 0.0, 0.0);
    CliffordCore::Vector3<double> b(0.0, 1.0, 0.0);
    auto dot = CliffordCore::dot_product(a, b);
    auto wedge = CliffordCore::wedge_product(a, b);
    auto product = CliffordCore::geometric_product(a, b);
    std::cout << dot.value << "\n";
    std::cout << wedge.xy << "\n";
    std::cout << product.scalar.value << "\n";
    return 0;
}
