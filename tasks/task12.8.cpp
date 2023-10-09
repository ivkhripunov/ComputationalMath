//
// Created by ivankhripunov on 09.10.23.
//

#include <gtest/gtest.h>
#include "NonLinearMPI.h"

/*
 * Замечание: для одной функции phi при разных начальных приближениях МПИ все равно сходится к одному из корней (а не к другому).
 * Чтобы найти другой корень, необходимо выразить х другим способом.
 */
TEST(TASK12_8, MPI) {
    const scalar sqr2 = std::sqrt(2);
    const auto equationFunc = [&](const scalar x) {
        return x * std::exp(-x * x) - 1 / (2 * std::sqrt(2 * M_E));
    };

    const Segment segment{0, 5};
    const indexType pointsCount = 100;

    const std::vector<scalar> localizedRoots = localizeRoots(equationFunc, pointsCount, segment);

    const scalar tolerance = 1e-3;
    const indexType maxIterationCount = 100;

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


/*
 * Наткнулся на проблему, что если взять маленький шаг, то выполняется критерий остановы и цикл завершается,
 * хотя фактической сходимости нет. Нужно доработать критерий остановы.
 */
TEST(TASK12_8, RELAXATION) {
    const auto equationFunc = [&](const scalar x) {
        return x * std::exp(-x * x) - 1 / (2 * std::sqrt(2 * M_E));
    };

    const Segment segment{0, 5};
    const indexType pointsCount = 100;

    const std::vector<scalar> localizedRoots = localizeRoots(equationFunc, pointsCount, segment);

    const scalar tolerance = 1e-3;
    const indexType maxIterationCount = 100;

    const scalar initialGuess1 = localizedRoots[0];
    const scalar initialGuess2 = localizedRoots[1];

    const scalar step = 1;
    const Result root1 = solveWithRelaxation(equationFunc, initialGuess1, step, tolerance, maxIterationCount);
    const Result root2 = solveWithRelaxation(equationFunc, initialGuess2, step, tolerance, maxIterationCount);
    const scalar referenceRoot1 = 0.226;
    const scalar referenceRoot2 = 1.359;

    ASSERT_TRUE(root1.converged);
    ASSERT_TRUE(root2.converged);
    ASSERT_NEAR(root1.value, referenceRoot1, tolerance);
    ASSERT_NEAR(root2.value, referenceRoot2, tolerance);
}
