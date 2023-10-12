//
// Created by ivankhripunov on 12.10.23.
//

#include <gtest/gtest.h>
#include "Chebyshev.h"

TEST(CHEBYSHEV, DEFAULT_SEGMENT1) {

    const indexType N = 0;
    std::array<scalar, N> chebRootsDefaultSegment = calcChebyshevRootsDefaultSegment<N>();
}

TEST(CHEBYSHEV, DEFAULT_SEGMENT2) {

    const scalar tolerance = 1e-13;
    const indexType N = 2;
    std::array<scalar, N> chebRootsDefaultSegment = calcChebyshevRootsDefaultSegment<N>();
    std::array<scalar, N> referenceRoots = {-1 / std::sqrt(2), 1 / std::sqrt(2)};

    for (indexType i = 0; i < N; ++i) {
        ASSERT_NEAR(chebRootsDefaultSegment[i], referenceRoots[i], tolerance);
    }
}

TEST(CHEBYSHEV, DEFAULT_SEGMENT3) {

    const scalar tolerance = 1e-13;
    const indexType N = 3;
    std::array<scalar, N> chebRootsDefaultSegment = calcChebyshevRootsDefaultSegment<N>();
    std::array<scalar, N> referenceRoots = {-std::sqrt(3) / 2, 0, std::sqrt(3) / 2};

    for (indexType i = 0; i < N; ++i) {
        ASSERT_NEAR(chebRootsDefaultSegment[i], referenceRoots[i], tolerance);
    }
}

TEST(CHEBYSHEV, DEFAULT_SEGMENT5) {

    const scalar tolerance = 1e-13;
    const indexType N = 5;
    std::array<scalar, N> chebRootsDefaultSegment = calcChebyshevRootsDefaultSegment<N>();
    std::array<scalar, N> referenceRoots = {-std::sqrt(2.5 + std::sqrt(5) / 2) / 2,
                                            -std::sqrt(2.5 - std::sqrt(5) / 2) / 2,
                                            0,
                                            std::sqrt(2.5 - std::sqrt(5) / 2) / 2,
                                            std::sqrt(2.5 + std::sqrt(5) / 2) / 2};

    for (indexType i = 0; i < N; ++i) {
        ASSERT_NEAR(chebRootsDefaultSegment[i], referenceRoots[i], tolerance);
    }
}

TEST(CHEBYSHEV, CUSTOM_SEGMENT0) {

    const scalar tolerance = 1e-13;
    const indexType N = 3;
    const Segment segment = {-1, 1};
    std::array<scalar, N> chebRootsCustomSegment = calcChebyshevRootsCustomSegment<N>(segment);
    std::array<scalar, N> referenceRoots = {-std::sqrt(3) / 2, 0, std::sqrt(3) / 2};

    for (indexType i = 0; i < N; ++i) {
        ASSERT_NEAR(chebRootsCustomSegment[i], referenceRoots[i], tolerance);
    }
}

TEST(CHEBYSHEV, CUSTOM_SEGMENT1) {

    const scalar tolerance = 1e-13;
    const indexType N = 3;
    const Segment segment = {0, 2};
    std::array<scalar, N> chebRootsCustomSegment = calcChebyshevRootsCustomSegment<N>(segment);
    std::array<scalar, N> referenceRoots = {-std::sqrt(3) / 2 + 1, 1, std::sqrt(3) / 2 + 1};

    for (indexType i = 0; i < N; ++i) {
        ASSERT_NEAR(chebRootsCustomSegment[i], referenceRoots[i], tolerance);
    }
}

TEST(CHEBYSHEV, CUSTOM_SEGMENT2) {

    const scalar tolerance = 1e-13;
    const indexType N = 5;
    const Segment segment = {0, 2};
    std::array<scalar, N> chebRootsCustomSegment = calcChebyshevRootsCustomSegment<N>(segment);
    std::array<scalar, N> referenceRoots = {-std::sqrt(2.5 + std::sqrt(5) / 2) / 2 + 1,
                                            -std::sqrt(2.5 - std::sqrt(5) / 2) / 2 + 1,
                                            1,
                                            std::sqrt(2.5 - std::sqrt(5) / 2) / 2 + 1,
                                            std::sqrt(2.5 + std::sqrt(5) / 2) / 2 + 1};

    for (indexType i = 0; i < N; ++i) {
        ASSERT_NEAR(chebRootsCustomSegment[i], referenceRoots[i], tolerance);
    }
}