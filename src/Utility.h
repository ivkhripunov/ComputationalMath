//
// Created by ivankhripunov on 11.10.23.
//

#ifndef COMPUTATIONALMATH_UTILITY_H
#define COMPUTATIONALMATH_UTILITY_H

#include <cmath>
#include <Eigen/Dense>

using scalar = double;
using indexType = std::size_t;
using Vector2d = Eigen::Vector2d;

template<indexType N>
using array = std::array<scalar, N>;

struct Segment {
    scalar begin;
    scalar end;
};

template <typename T, indexType N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr) {
    os << "[ ";
    for (const auto& element : arr) {
        os << element << " ";
    }
    os << "]";
    return os;
}

#endif //COMPUTATIONALMATH_UTILITY_H
