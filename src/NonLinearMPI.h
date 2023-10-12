//
// Created by ivankhripunov on 09.10.23.
//

#ifndef COMPUTATIONALMATH_NONLINEARMPI_H
#define COMPUTATIONALMATH_NONLINEARMPI_H

#include <vector>
#include "Utility.h"

struct Result {
    scalar value;
    bool converged;
};

struct ResultVector {
    Vector2d value;
    bool converged;
};

[[nodiscard]]
std::vector<scalar> inline calcGrid(const indexType pointsCount, const Segment &segment) {

    const auto gridSegmentCount = static_cast<scalar>(pointsCount - 1);
    const scalar step = (segment.end - segment.begin) / gridSegmentCount;

    std::vector<scalar> grid(pointsCount);

    for (indexType i = 0; i < pointsCount; ++i) {
        grid[i] = segment.begin + static_cast<scalar>(i) * step;
    }

    return grid;
}

//TODO: add localization with tolerance
[[nodiscard]]
std::vector<scalar>
localizeRoots(const auto &equationFunc, const Segment &segment, const indexType pointsCount) noexcept {

    std::vector<scalar> result;
    const std::vector<scalar> grid = calcGrid(pointsCount, segment);
    scalar valueInPointOld = equationFunc(grid[0]);

    for (indexType i = 0; i < pointsCount - 1; ++i) {

        const scalar valueInPointNew = equationFunc(grid[i + 1]);

        if (valueInPointOld * valueInPointNew <= 0) {
            result.push_back((grid[i] + grid[i + 1]) / 2);
        }

        valueInPointOld = valueInPointNew;
    }

    return result;
}

[[nodiscard]]
scalar inline
dichotomySolver(const auto &equationFunc, const Segment &initialSegment,
                const scalar tolerance) noexcept {

    scalar midPoint;
    Segment iterationSegment = initialSegment;

    while (iterationSegment.end - iterationSegment.begin >= tolerance) {

        midPoint = (iterationSegment.begin + iterationSegment.end) / 2;

        if (equationFunc(midPoint) * equationFunc(iterationSegment.begin) < 0) {
            iterationSegment.end = midPoint;
        } else { iterationSegment.begin = midPoint; }

    }

    return midPoint;
}


//Обязательно подставляем уже итеративную функцию!!!
[[nodiscard]]
Result inline solveWithMPI(const auto &iterativeFunc, const scalar initialGuess, const scalar tolerance,
                           const indexType maxIterationCount) noexcept {

    scalar x_new = initialGuess;
    indexType counter = 0;
    bool continueCriteria = true;

    while (continueCriteria) {
        const scalar x_old = x_new;
        x_new = iterativeFunc(x_old);

        counter++;
        continueCriteria = (counter < maxIterationCount) && (std::abs(x_new - x_old) > tolerance);
    }

    return {x_new, counter < maxIterationCount};
}

//Обязательно подставляем уже итеративную функцию!!!
[[nodiscard]]
ResultVector inline solveVectorWithMPI(const auto &iterativeFunc, Vector2d initialGuess,
                                       const scalar tolerance,
                                       const indexType maxIterationCount) noexcept {

    Vector2d x_new = initialGuess;
    indexType counter = 0;
    bool continueCriteria = true;

    while (continueCriteria) {
        const Vector2d x_old = x_new;
        x_new = iterativeFunc(x_old);

        counter++;
        continueCriteria = (counter < maxIterationCount) && ((x_new - x_old).norm() > tolerance);
    }

    return {x_new, counter < maxIterationCount};
}

template<bool positiveDerivative>
[[nodiscard]]
Result inline
solveWithRelaxation(const auto &equationFunc, const scalar initialGuess, const scalar step,
                    const indexType maxIterationCount, const scalar tolerance) noexcept {

    const scalar sign = positiveDerivative ? -1 : 1;
    scalar x_new = initialGuess;
    indexType counter = 0;

    const auto iterativeFunc = [&](const scalar x) { return sign * step * equationFunc(x) + x; };

    while (counter < maxIterationCount) {
        const scalar x_old = x_new;
        x_new = iterativeFunc(x_old);

        counter++;
    }

    return {x_new, std::abs(equationFunc(x_new)) < tolerance};
}

[[nodiscard]]
Result inline
solveWithNewton(const auto &equationFunc, const auto &derivativeFunc, const scalar initialGuess,
                const indexType maxIterationCount, const scalar tolerance) noexcept {

    scalar x_new = initialGuess;
    indexType counter = 0;
    const auto iterativeFunc = [&](const scalar x) { return x - equationFunc(x) / derivativeFunc(x); };


    while (counter < maxIterationCount) {
        const scalar x_old = x_new;
        x_new = iterativeFunc(x_old);

        counter++;
    }

    return {x_new, std::abs(equationFunc(x_new)) < tolerance};
}

[[nodiscard]]
ResultVector inline solveVectorWithNewton(const auto &equationFunc, const auto &Jacobi, Vector2d initialGuess,
                                          const scalar tolerance,
                                          const indexType maxIterationCount) noexcept {

    Vector2d x_new = initialGuess;
    indexType counter = 0;
    bool continueCriteria = true;

    while (continueCriteria) {
        const Vector2d x_old = x_new;
        const Eigen::Matrix<scalar, 2, 2> jacobi = Jacobi(x_old);
        const Eigen::Matrix<scalar, 2, 2> jacobiInv = jacobi.inverse();
        x_new = x_old - jacobiInv * equationFunc(x_old);

        counter++;
        continueCriteria = (counter < maxIterationCount) && ((x_new - x_old).norm() > tolerance);
    }

    return {x_new, counter < maxIterationCount};
}

#endif //COMPUTATIONALMATH_NONLINEARMPI_H
