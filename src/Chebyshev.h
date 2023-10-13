//
// Created by ivankhripunov on 12.10.23.
//

#ifndef COMPUTATIONALMATH_CHEBYSHEV_H
#define COMPUTATIONALMATH_CHEBYSHEV_H

#include <algorithm>
#include "Utility.h"

template<indexType N>
[[nodiscard]]
constexpr std::array<scalar, N> calcChebyshevRootsDefaultSegment() noexcept {

    const auto &rootFunc = [](const indexType i) {
        return std::cos((2 * static_cast<scalar>(i) + 1) * M_PI / (2 * N));
    };

    std::array<scalar, N> result;

    for (indexType i = 0; i < N; ++i) {
        result[i] = rootFunc(i);
    }

    std::sort(result.begin(), result.end());

    return result;
}

template<indexType N>
[[nodiscard]]
constexpr
std::array<scalar, N>
transformToCustomSegment(const std::array<scalar, N> &rootsDefaultSegment,
                                     const Segment &segment) noexcept {

    const scalar aPlusBhalf = (segment.begin + segment.end) / 2;
    const scalar aMinusBhalf = (segment.end - segment.begin) / 2;

    const auto &transformRoot = [&](const scalar root) {
        return aPlusBhalf + aMinusBhalf * root;
    };

    std::array<scalar, N> result;

    for (indexType i = 0; i < N; ++i) {
        result[i] = transformRoot(rootsDefaultSegment[i]);
    }

    return result;
}

template<indexType N>
[[nodiscard]]
constexpr
std::array<scalar, N>
calcChebyshevRootsCustomSegment(const Segment &segment) noexcept {

    return transformToCustomSegment(calcChebyshevRootsDefaultSegment<N>(), segment);
}

#endif //COMPUTATIONALMATH_CHEBYSHEV_H
