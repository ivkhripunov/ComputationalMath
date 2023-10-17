//
// Created by ivankhripunov on 16.10.23.
//

#include <algorithm>
#include "NonLinearMPI.h"
#include <gtest/gtest.h>
#include "Spline.h"

TEST(SPLINE, SET1) {
    const scalar tolerance = 1e-11;
    const indexType N = 5;
    const Segment segment = {0, 10};
    const auto x = calcGrid(N, segment);

    const auto function = [](const scalar x) { return std::exp(x); };
    std::vector<scalar> y(N);

    for (indexType i = 0; i < N; ++i) {
        y[i] = function(x[i]);
    }

    CubicSpline<scalar, scalar> splineInterpolator(x, y);

    for (indexType i = 0; i < 1; ++i) {
        ASSERT_NEAR(splineInterpolator.interpolate(x[i]), y[i], tolerance);
    }

    const indexType checkerN = 1000;
    const auto checkerGrid = calcGrid(checkerN, segment);

    scalar errorMax = 0;
    for (const auto &point: checkerGrid) {

        const scalar interpolationResult = splineInterpolator.interpolate(point);
        const scalar referenceResult  = function(point);
        const scalar relativeError = std::abs((interpolationResult - referenceResult));

        if (errorMax < relativeError) {
            errorMax = relativeError;
            //std::cout << point << " ";
        }
    }

    std::cout << errorMax << std::endl;
}

