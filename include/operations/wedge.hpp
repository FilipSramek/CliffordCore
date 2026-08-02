#pragma once

#include <type_traits>
#include "../vector3.hpp"
#include "../bivector3.hpp"

namespace CliffordCore
{

template<typename T>
constexpr Bivector3<T> operator^(const Vector3<T>& a, const Vector3<T>& b) {
    return Bivector3<T>(
        a.x * b.y - a.y * b.x,  // xy component
        a.x * b.z - a.z * b.x,  // xz component
        a.y * b.z - a.z * b.y   // yz component
    );
};
}