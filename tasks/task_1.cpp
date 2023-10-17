//
// Created by ivankhripunov on 11.10.23.
//

#include <gtest/gtest.h>
#include "NonLinearMPI.h"
#include <Eigen/Dense>

//стандартный MPI разваливается!!! Все плохо в 0 из-за корня
TEST(TASK1, MPI) {

    const auto iterativeFunc1 = [](const Vector2d &vector) {
        return Vector2d(
                std::sqrt(1 - vector.y() * vector.y()),
                std::tan(vector.x())
        );
    };


    const auto iterativeFunc2 = [](const Vector2d &vector) {
        return Vector2d(
                -std::sqrt(1 - vector.y() * vector.y()),
                std::tan(vector.x())
        );
    };

//    const Vector2d initialGuess1 = {0, 0};
//    const Vector2d initialGuess2 = {0, 0};

    const Vector2d initialGuess1 = {M_SQRT2 / 2, M_SQRT2 / 2};
    const Vector2d initialGuess2 = {-M_SQRT2 / 2, -M_SQRT2 / 2};

    const scalar tolerance = 1e-6;
    const indexType maxIterationCount = 20;

    const ResultVector root1 = solveVectorWithMPI(iterativeFunc1, initialGuess1, tolerance, maxIterationCount);
    const ResultVector root2 = solveVectorWithMPI(iterativeFunc2, initialGuess2, tolerance, maxIterationCount);

    const Vector2d referenceRoot1 = {0.649889, 0.760029};
    const Vector2d referenceRoot2 = -referenceRoot1;

    ASSERT_TRUE(root1.converged);
    ASSERT_TRUE(root2.converged);
    ASSERT_LE((root1.value - referenceRoot1).norm(), tolerance);
    ASSERT_LE((root2.value - referenceRoot2).norm(), tolerance);
}

TEST(TASK1, RELAXATION) {

    const auto equationFunc = [](const Vector2d &vector) {
        return Vector2d(
                vector.squaredNorm() - 1,
                vector.y() - std::tan(vector.x())
        );
    };


    const Vector2d initialGuess1 = {M_SQRT2 / 2, M_SQRT2 / 2};
    const Vector2d initialGuess2 = {-M_SQRT2 / 2, -M_SQRT2 / 2};

    const scalar tolerance = 1e-6;
    const indexType maxIterationCount = 100;

    const scalar step1 = -0.2;
    const scalar step2 = 0.2;

    const ResultVector root1 = solveVectorWithRelaxation(equationFunc, initialGuess1, step1, maxIterationCount, tolerance);
    const ResultVector root2 = solveVectorWithRelaxation(equationFunc, initialGuess2, step2, maxIterationCount, tolerance);

    const Vector2d referenceRoot1 = {0.649889, 0.760029};
    const Vector2d referenceRoot2 = -referenceRoot1;

    ASSERT_TRUE(root1.converged);
    ASSERT_TRUE(root2.converged);
    ASSERT_LE((root1.value - referenceRoot1).norm(), tolerance);
    ASSERT_LE((root2.value - referenceRoot2).norm(), tolerance);
}

TEST(TASK1, NEWTON) {

    const auto equationFunc = [](const Vector2d &vector) {
        return Vector2d(
                vector.squaredNorm() - 1,
                vector.y() - std::tan(vector.x())
        );
    };


    const auto jacobiFunc = [](const Vector2d &vector) {
        const scalar cosx = std::cos(vector.x());
        Eigen::Matrix<scalar, 2, 2> jacobi;
        jacobi(0, 0) = 2 * vector.x();
        jacobi(0, 1) = 2 * vector.y();
        jacobi(1, 0) = -1 / (cosx * cosx);
        jacobi(1, 1) = 1;

        return jacobi;
    };

//    const Vector2d initialGuess1 = {1, 1};
//    const Vector2d initialGuess2 = {-1, -1};

    const Vector2d initialGuess1 = {M_SQRT2 / 2, M_SQRT2 / 2};
    const Vector2d initialGuess2 = {-M_SQRT2 / 2, -M_SQRT2 / 2};

    const scalar tolerance = 1e-6;
    const indexType maxIterationCount = 5;

    const ResultVector root1 = solveVectorWithNewton(equationFunc, jacobiFunc, initialGuess1, tolerance, maxIterationCount);
    const ResultVector root2 = solveVectorWithNewton(equationFunc, jacobiFunc, initialGuess2, tolerance, maxIterationCount);

    const Vector2d referenceRoot1 = {0.649889, 0.760029};
    const Vector2d referenceRoot2 = -referenceRoot1;

    ASSERT_TRUE(root1.converged);
    ASSERT_TRUE(root2.converged);
    ASSERT_LE((root1.value - referenceRoot1).norm(), tolerance);
    ASSERT_LE((root2.value - referenceRoot2).norm(), tolerance);
}