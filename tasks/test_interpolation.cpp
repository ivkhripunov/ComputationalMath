//
// Created by ivankhripunov on 12.10.23.
//

#include <gtest/gtest.h>
#include "Interpolation.h"

TEST(interpolationTest, Test_1) {
    const std::size_t N = 10;
    srand((unsigned) time(NULL));
    std::array<scalar, N> xArr, yArr;
    xArr[0] = static_cast<scalar>(std::rand()) / RAND_MAX * 4;
    yArr[0] = static_cast<scalar>(std::rand()) / RAND_MAX * 40;
    for (int i = 1; i < N; i++) {
        xArr[i] = xArr[i - 1] + static_cast<scalar>(std::rand()) / RAND_MAX * 4;
        yArr[i] = yArr[i - 1] + static_cast<scalar>(std::rand()) / RAND_MAX * 40;
    }
    for (int i = 0; i < N; i++) {
        std::cout << xArr[i] << " " << yArr[i] << std::endl;
    }
    NewtonInterpolator<scalar, scalar, N> interpolator(xArr, yArr);
    for (size_t i = 0; i < N; ++i) {
        ASSERT_NEAR((yArr[i] - interpolator.interpolate(xArr[i])) / yArr[i], 0, 1e-12);
    }
}
