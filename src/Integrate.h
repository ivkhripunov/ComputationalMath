//
// Created by ivankhripunov on 17.10.23.
//

#ifndef COMPUTATIONALMATH_INTEGRATE_H
#define COMPUTATIONALMATH_INTEGRATE_H

#include "Utility.h"
#include "Chebyshev.h"

template<typename T>
[[nodiscard]] T integrateSegmentTrapecia(const auto &function, const Segment &segment) noexcept {

    return (function(segment.begin) + function(segment.end)) * (segment.end - segment.begin) / 2;
}

template<typename T>
[[nodiscard]] T
integrateGridTrapecia(const auto &function, const std::vector<T> &grid) noexcept {


    scalar result = 0;
    for (indexType i = 0; i < grid.size() - 1; ++i) {
        result += integrateSegmentTrapecia<scalar>(function, {grid[i], grid[i + 1]});
    }
    return result;
}

template<typename T>
[[nodiscard]] T integrateSegmentSimpson(const auto &function, const Segment &segment) noexcept {

    return (function(segment.begin) + 4 * function((segment.begin + segment.end) / 2) + function(segment.end))
           * (segment.end - segment.begin) / 6;
}

template<typename T>
[[nodiscard]] T
integrateGridSimpson(const auto &function, const std::vector<T> &grid) noexcept {

    scalar result = 0;
    for (indexType i = 0; i < grid.size() - 1; ++i) {
        result += integrateSegmentSimpson<scalar>(function, {grid[i], grid[i + 1]});
    }
    return result;
}


[[nodiscard]] scalar
integrateSegmentGauss(const auto &function, const Segment segment) noexcept {

    const scalar aPlusBhalf = (segment.begin + segment.end) / 2;
    const scalar aMinusBhalf = (segment.end - segment.begin) / 2;
    const auto &transformPoint = [&](const scalar root) {

        return aPlusBhalf + aMinusBhalf * root;
    };

    const indexType N = 3;
    const scalar sqr = std::sqrt(0.6);
    const std::array<scalar, N> points = transformToCustomSegment<3>({-sqr, 0, sqr}, segment);
    const std::array<scalar, N> weights = {5. / 9, 8. / 9, 5. / 9};

    scalar result = 0;
    for (indexType i = 0; i < N; ++i) {
        result += function(points[i]) * weights[i] * aMinusBhalf;
    }

    return result;
}


#endif //COMPUTATIONALMATH_INTEGRATE_H
