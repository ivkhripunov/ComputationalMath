//
// Created by ivankhripunov on 17.10.23.
//

#include <algorithm>
#include "NonLinearMPI.h"
#include <gtest/gtest.h>
#include "Spline.h"
#include "Interpolation.h"

TEST(TASK9_32, EXTRAPOLATION) {

    const indexType N = 10;
    const std::array<scalar, N> yearsArray = {1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000};
    const std::array<scalar, N> populationArray = {92228496, 106021537, 123202624, 132164569, 151325798, 179323175, 203211926,
                                                   226545805, 248709873, 281421906};

    const std::vector<scalar> yearsVector = {1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000};
    const std::vector<scalar> populationVector = {92228496, 106021537, 123202624, 132164569, 151325798, 179323175, 203211926,
                                              226545805, 248709873, 281421906};


    const NewtonInterpolator<scalar, scalar, 10> newtonInterpolator(yearsArray, populationArray);
    const CubicSpline<scalar, scalar> splineInterpolator(yearsVector, populationVector);

    const scalar referenceResult = 308745538;

    std::cout << splineInterpolator.extrapolateRight(2010) << std::endl;

    std::cout << newtonInterpolator.interpolate(2010);
}
