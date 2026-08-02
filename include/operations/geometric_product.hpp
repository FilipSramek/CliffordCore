#pragma once

#include <type_traits>
#include "../vector3.hpp"
#include "../multivector3.hpp"
#include "dot_product.hpp"
#include "wedge_product.hpp"

namespace CliffordCore
{
    template<typename T>
    /**
     * @brief Computes the geometric product of two 3D vectors, resulting in a multivector.
     * @param a The first vector.
     * @param b The second vector.
     * @return The resulting multivector from the geometric product = dot_product(a, b) + wedge_product(a, b) so that the scalar part is the dot product and the bivector part is the wedge product.
     */
    constexpr Multivector3<T> geometric_product(const Vector3<T>& a, const Vector3<T>& b) {
        return Multivector3<T>(
            dot_product(a, b),
            vector3<T>(0, 0, 0),
            wedge_product(a, b),
            trivector3<T>(0)
        );
    }
} // namespace CliffordCore