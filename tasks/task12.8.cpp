//
// Created by ivankhripunov on 09.10.23.
//

#include <gtest/gtest.h>
#include "NonLinearMPI.h"
/*
 * Замечание: для одной функции phi при разных начальных приближениях МПИ все равно сходится к одному из корней (а не к другому).
 * Чтобы найти другой корень, необходимо выразить х другим способом. Это также свзяано с условием сходимости.
 */
TEST(TASK12_8, MPI) {
    const scalar sqr2 = std::sqrt(2);
    const auto equationFunc = [&](const scalar x) {
        return x * std::exp(-x * x) - 1 / (2 * std::sqrt(2 * M_E));
    };

    const Segment segment{0, 5};
    const indexType pointsCount = 10;

    const std::vector<scalar> localizedRoots = localizeRoots(equationFunc, segment, pointsCount);

    const scalar tolerance = 1e-3;
    const indexType maxIterationCount = 5;

    const auto iterativeFunc1 = [&](const scalar x) {
        return std::exp(x * x - 0.5) / (2 * sqr2);
    };

    const auto iterativeFunc2 = [&](const scalar x) {
        return std::sqrt(std::log(2 * x * std::sqrt(2 * M_E)));
    };

    const scalar initialGuess1 = localizedRoots[0];
    const scalar initialGuess2 = localizedRoots[1];

    const Result root1 = solveWithMPI(iterativeFunc1, initialGuess1, tolerance, maxIterationCount);
    const Result root2 = solveWithMPI(iterativeFunc2, initialGuess2, tolerance, maxIterationCount);
    const scalar referenceRoot1 = 0.226;
    const scalar referenceRoot2 = 1.359;

    ASSERT_TRUE(root1.converged);
    ASSERT_TRUE(root2.converged);
    ASSERT_NEAR(root1.value, referenceRoot1, tolerance);
    ASSERT_NEAR(root2.value, referenceRoot2, tolerance);
}

TEST(TASK12_8, RELAXATION) {
    const auto equationFunc = [](const scalar x) {
        return x * std::exp(-x * x) - 1 / (2 * std::sqrt(2 * M_E));
    };

    const Segment segment{0, 5};
    const indexType pointsCount = 10;

    const std::vector<scalar> localizedRoots = localizeRoots(equationFunc, segment, pointsCount);

    const scalar tolerance = 1e-3;
    const indexType maxIterationCount = 12;
    const scalar step1 = 2;
    const scalar step2 = 2;

    const scalar initialGuess1 = localizedRoots[0];
    const scalar initialGuess2 = localizedRoots[1];

    const Result root1 = solveWithRelaxation<true>(equationFunc, initialGuess1, step1, maxIterationCount, tolerance);
    const Result root2 = solveWithRelaxation<false>(equationFunc, initialGuess2, step2, maxIterationCount, tolerance);
    const scalar referenceRoot1 = 0.226;
    const scalar referenceRoot2 = 1.359;

    ASSERT_TRUE(root1.converged);
    ASSERT_TRUE(root2.converged);
    ASSERT_NEAR(root1.value, referenceRoot1, tolerance);
    ASSERT_NEAR(root2.value, referenceRoot2, tolerance);
}

TEST(TASK12_8, DICHOTOMY) {
    const auto equationFunc = [](const scalar x) {
        return x * std::exp(-x * x) - 1 / (2 * std::sqrt(2 * M_E));
    };
    const scalar max = 1 / std::sqrt(2 * M_E);
    const Segment segment1{0, max};
    const Segment segment2{max, 5};
    const scalar tolerance = 1e-3;

    const scalar root1 = dichotomySolver(equationFunc, segment1, tolerance);
    const scalar root2 = dichotomySolver(equationFunc, segment2, tolerance);
    const scalar referenceRoot1 = 0.226;
    const scalar referenceRoot2 = 1.359;

    ASSERT_NEAR(root1, referenceRoot1, tolerance);
    ASSERT_NEAR(root2, referenceRoot2, tolerance);
}

TEST(TASK12_8, NEWTON) {
    const auto equationFunc = [](const scalar x) {
        return x * std::exp(-x * x) - 1 / (2 * std::sqrt(2 * M_E));
    };

    const auto derivativeFunc = [](const scalar x) {
        const scalar xSqr = x * x;
        return std::exp(-xSqr) * (1 - 2 * xSqr);
    };

    const Segment segment{0, 5};
    const indexType pointsCount = 7;

    const std::vector<scalar> localizedRoots = localizeRoots(equationFunc, segment, pointsCount);

    const scalar tolerance = 1e-3;
    const indexType maxIterationCount = 5;

    const scalar initialGuess1 = localizedRoots[0];
    const scalar initialGuess2 = localizedRoots[1];

    const Result root1 = solveWithNewton(equationFunc, derivativeFunc, initialGuess1, maxIterationCount, tolerance);
    const Result root2 = solveWithNewton(equationFunc, derivativeFunc, initialGuess2, maxIterationCount, tolerance);
    const scalar referenceRoot1 = 0.226;
    const scalar referenceRoot2 = 1.359;

    ASSERT_TRUE(root1.converged);
    ASSERT_TRUE(root2.converged);
    ASSERT_NEAR(root1.value, referenceRoot1, tolerance);
    ASSERT_NEAR(root2.value, referenceRoot2, tolerance);
}
