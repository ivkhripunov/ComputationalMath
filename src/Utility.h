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

[[nodiscard]]
std::vector<scalar> inline calcGrid(const indexType pointsCount, const Segment &segment) {

    const auto gridSegmentCount = static_cast<scalar>(pointsCount - 1);
    const scalar step = (segment.end - segment.begin) / gridSegmentCount;

    std::vector<scalar> grid(pointsCount);

    for (indexType i = 0; i < pointsCount; ++i) {
        grid[i] = segment.begin + static_cast<scalar>(i) * step;
    }

    return grid;
}

template <typename T, indexType N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr) {
    os << "[ ";
    for (const auto& element : arr) {
        os << element << " ";
    }
    os << "]";
    return os;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    os << "[ ";
    for (const auto& element : vec) {
        os << element << " ";
    }
    os << "]";
    return os;
}

#endif //COMPUTATIONALMATH_UTILITY_H
