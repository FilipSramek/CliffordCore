#pragma once

#include "../scalar.hpp"
#include "../vector3.hpp"

namespace CliffordCore
{

template<typename T>
constexpr Scalar<T> operator*(const Vector3<T>& a, const Vector3<T>& b) {
    return Scalar<T>(a.x * b.x + a.y * b.y + a.z * b.z);
}

} // namespace CliffordCore