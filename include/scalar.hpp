#pragma once

#include <type_traits>

namespace CliffordCore {

template <typename T>

class Scalar {    
    static_assert(std::is_arithmetic<T>::value, "Scalar can only be instantiated with numeric types.");

public:
    T value;

    constexpr Scalar() : value(0) {}                            // Default constructor initializes value to zero
    constexpr Scalar(T val) : value(val) {}                     // Constructor initializes value to the provided argument

    constexpr Scalar operator+(const Scalar& other) const {     // Addition operator overload
        return Scalar(value + other.value);
    }

    constexpr Scalar operator-(const Scalar& other) const {     // Subtraction operator overload
        return Scalar(value - other.value);
    }

    constexpr Scalar operator*(const Scalar& other) const {     // Multiplication operator overload
        return Scalar(value * other.value);
    }

    constexpr Scalar operator/(const Scalar& other) const {     // Division operator overload
        return Scalar(value / other.value);                     
    }

    constexpr Scalar operator-() const {                        // Unary negation operator overload 
        return Scalar(-value);
    }
};
}