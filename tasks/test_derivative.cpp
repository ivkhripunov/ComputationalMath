//
// Created by ivankhripunov on 11.10.23.
//

#include <gtest/gtest.h>
#include "Derivative.h"

TEST(differentiationtest_test, Test_1)
{
    const std::size_t N = 2;
    std::array<double, 2> hCoeff{-1, 1};
    const double trueCentralCoeff = 0.0;
    const std::array<double, 2> trueOtherCoeff {-0.5, 0.5};
    calcDerivativeCoeff<double, N> test_answer = calcDerivativeCoef<double, N, 1>(hCoeff);
    for (size_t i = 0; i < N; ++i)
    {
        ASSERT_NEAR(trueOtherCoeff[i], test_answer.otherPointsCoeff[i], 1e-13);
    }
    ASSERT_NEAR(trueCentralCoeff, test_answer.centralPointCoeff, 1e-13);
}

TEST(differentiationtest_test, Test_2)
{
    const std::size_t N = 2;
    std::array<double, 2> hCoeff{1, 2};
    const double trueCentralCoeff = -1.5;
    const std::array<double, 2> trueOtherCoeff {2.0, -0.5};
    calcDerivativeCoeff<double, N> test_answer = calcDerivativeCoef<double, N, 1>(hCoeff);
    for (size_t i = 0; i < N; ++i)
    {
        ASSERT_NEAR(trueOtherCoeff[i], test_answer.otherPointsCoeff[i], 1e-13);
    }
    ASSERT_NEAR(trueCentralCoeff, test_answer.centralPointCoeff, 1e-13);
}

TEST(differentiationtest_test, Test_3)
{
    const std::size_t N = 2;
    std::array<double, 2> hCoeff{-1, 1};
    const double trueCentralCoeff = -2;
    const std::array<double, 2> trueOtherCoeff {1, 1};
    calcDerivativeCoeff<double, N> test_answer = calcDerivativeCoef<double, N, 2>(hCoeff);
    for (size_t i = 0; i < N; ++i)
    {
        ASSERT_NEAR(trueOtherCoeff[i], test_answer.otherPointsCoeff[i], 1e-13);
    }
    ASSERT_NEAR(trueCentralCoeff, test_answer.centralPointCoeff, 1e-13);
}