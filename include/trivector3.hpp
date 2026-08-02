#pragma once

#include <type_traits>

namespace CliffordCore
{

template<typename T>
class Trivector3
{
    static_assert(std::is_arithmetic<T>::value, "Trivector3 can only be instantiated with numeric types.");

public:
    T e123;
};

} // namespace CliffordCore