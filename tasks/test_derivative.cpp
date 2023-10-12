//
// Created by ivankhripunov on 11.10.23.
//

#include <gtest/gtest.h>
#include "Derivative.h"

TEST(DERIVATIVE, SET1) {
    const indexType N = 2;
    const std::array<scalar, 2> hCoeff{-1, 1};
    const scalar referenceCentralCoeff = 0;
    const std::array<scalar, 2> trueOtherCoeff{-0.5, 0.5};
    const calcDerivativeCoeff<scalar, N> test_answer = calcDerivativeCoef<scalar, N, 1>(hCoeff);
    for (size_t i = 0; i < N; ++i) {
        ASSERT_NEAR(trueOtherCoeff[i], test_answer.otherPointsCoeff[i], 1e-13);
    }
    ASSERT_NEAR(referenceCentralCoeff, test_answer.centralPointCoeff, 1e-13);
}

TEST(DERIVATIVE, SET2) {
    const indexType N = 2;
    const std::array<scalar, 2> hCoeff{1, 2};
    const scalar referenceCentralCoeff = -1.5;
    const std::array<scalar, 2> trueOtherCoeff{2.0, -0.5};
    const calcDerivativeCoeff<scalar, N> test_answer = calcDerivativeCoef<scalar, N, 1>(hCoeff);
    for (size_t i = 0; i < N; ++i) {
        ASSERT_NEAR(trueOtherCoeff[i], test_answer.otherPointsCoeff[i], 1e-13);
    }
    ASSERT_NEAR(referenceCentralCoeff, test_answer.centralPointCoeff, 1e-13);
}

TEST(DERIVATIVE, SET3) {
    const indexType N = 2;
    const std::array<scalar, 2> hCoeff{-1, 1};
    const scalar referenceCentralCoeff = -2;
    const std::array<scalar, 2> trueOtherCoeff{1, 1};
    const calcDerivativeCoeff<scalar, N> test_answer = calcDerivativeCoef<scalar, N, 2>(hCoeff);
    for (size_t i = 0; i < N; ++i) {
        ASSERT_NEAR(trueOtherCoeff[i], test_answer.otherPointsCoeff[i], 1e-13);
    }
    ASSERT_NEAR(referenceCentralCoeff, test_answer.centralPointCoeff, 1e-13);
}