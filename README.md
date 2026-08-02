# CliffordCore

CliffordCore is a small, header-only C++ library for working with elements of 3D Clifford algebra. It provides basic types for scalars, vectors, bivectors, trivectors, multivectors, and rotors, along with common operations such as dot products, wedge products, geometric products, norms, inverses, duality, and reversals.

## Features

- 3D vector, bivector, trivector, multivector, and rotor types
- Arithmetic operators for the core algebraic objects
- Support for:
  - dot product
  - wedge product
  - geometric product
  - norm and squared norm
  - inverse
  - reverse
  - dual
- Header-only design for straightforward integration into C++ projects

## Project structure

- include/ — core headers and operation headers
- examples/ — example usage files
- tests/ — test cases and validation code
- docs/ — documentation and notes

## Usage

Include the headers you need from the include directory.

```cpp
#include <iostream>
#include "include/vector3.hpp"
#include "include/operations/dot_product.hpp"
#include "include/operations/wedge_product.hpp"
#include "include/operations/geometric_product.hpp"

int main() {
    CliffordCore::Vector3<double> a(1.0, 0.0, 0.0);
    CliffordCore::Vector3<double> b(0.0, 1.0, 0.0);

    auto dot = CliffordCore::dot_product(a, b);
    auto wedge = CliffordCore::wedge_product(a, b);
    auto product = CliffordCore::geometric_product(a, b);

    std::cout << "dot: " << dot.value << "\n";
    std::cout << "wedge: " << wedge.xy << ", " << wedge.xz << ", " << wedge.yz << "\n";
    return 0;
}
```

## Build and usage notes

This project is currently implemented as a set of header files, so it can be used directly by adding the repository's include directory to your compiler include path.

## Status

CliffordCore is an educational and experimental implementation of basic 3D Clifford algebra concepts. It is suitable for learning, prototyping, and exploring geometric algebra in C++.
