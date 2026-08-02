#pragma once

#include <type_traits>
#include "../vector3.hpp"
#include "dot_product.hpp"
#include "norm.hpp"

namespace CliffordCore
{
    template<typename T>
    /**
     * @brief Computes the inverse of a 3D vector.
     * @param v The vector for which to compute the inverse.
     * @return The inverse vector to the input vector.
     */
    constexpr Vector3<T> inverse(const Vector3<T>& v) {
        return v/squared_norm(v);
    }
} // namespace CliffordCore