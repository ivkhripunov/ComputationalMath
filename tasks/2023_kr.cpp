//
// Created by ivankhripunov on 12.10.23.
//

#include <gtest/gtest.h>
#include "Chebyshev.h"
#include "Derivative.h"
#include "Interpolation.h"
#include "NonLinearMPI.h"

//Обратная интерполяция
TEST(KR, TASK5) {

    const indexType N = 3;
    const Segment segment = {-2, -1};
    const array<N> xPoints = calcChebyshevRootsCustomSegment<N>(segment);
    const array<N> yPoints = {-0.5, -0.1, 0.3};

    const NewtonInterpolator<scalar, scalar, N> interpolator(yPoints, xPoints);

    std::cout << "Ответ:" << interpolator.interpolate(0) << std::endl;

    const NewtonInterpolator<scalar, scalar, N> checker(xPoints, yPoints);

    std::cout << checker.interpolate(interpolator.interpolate(0));
}

//Дифференцирование
TEST(KR, TASK4) {

    const indexType N = 2;
    const array<N> points = {-2, 1};
    const derivativeCoeff<scalar, N> coeff = calcDerivativeCoef<scalar, N, 1>(points);

    std::cout << "Side:   " << coeff.otherPointsCoeff << std::endl << "Central:   " << coeff.centralPointCoeff;
}

TEST(KR, TASK2) {

    const indexType N = 3;
    const array<N> xPoints = {-1, -0.3, 1};
    const array<N> yPoints = {std::sin(-1), std::sin(-0.3), std::sin(1)};

    const NewtonInterpolator<scalar, scalar, N> interpolator(xPoints, yPoints);


    std::cout << interpolator.interpolate(1) << std::endl;
    std::cout << interpolator.interpolate(3 / 2 + 0.3 * (-0.5)) << std::endl;
    std::cout << interpolator.interpolate(2) << std::endl;
}

TEST(KR, TASK7) {

    const indexType N = 2;
    const array<N> xPoints = {1, 2};
    const array<N> yPoints = {1, 0.5};

    const NewtonInterpolator<scalar, scalar, N> interpolator(xPoints, yPoints);
}