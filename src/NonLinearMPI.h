//
// Created by ivankhripunov on 09.10.23.
//

#ifndef COMPUTATIONALMATH_NONLINEARMPI_H
#define COMPUTATIONALMATH_NONLINEARMPI_H

#include <cmath>
#include <vector>

using scalar = double;
using indexType = std::size_t;

struct Segment {
    scalar begin;
    scalar end;
};

struct Result {
    scalar value;
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
std::vector<scalar> inline
localizeRoots(const auto &func, const indexType pointsCount, const Segment &segment) noexcept {

    std::vector<scalar> result;
    const std::vector<scalar> grid = calcGrid(pointsCount, segment);
    scalar valueInPointOld = func(grid[0]);

    for (indexType i = 0; i < pointsCount - 1; ++i) {

        const scalar valueInPointNew = func(grid[i + 1]);

        if (valueInPointOld * valueInPointNew <= 0) {
            result.push_back((grid[i] + grid[i + 1]) / 2);
        }

        valueInPointOld = valueInPointNew;
    }

    return result;
}

//RHS = 0
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

//TODO: fix relaxation
[[nodiscard]]
Result inline
solveWithRelaxation(const auto &equationFunc, const scalar initialGuess, const scalar step, const scalar tolerance,
                    const indexType maxIterationCount) noexcept {

    scalar x_new = initialGuess;
    indexType counter = 0;
    bool continueCriteria = true;

    const auto iterativeFunc = [&](const scalar x) { return step * equationFunc(x) + x; };

    while (continueCriteria) {
        const scalar x_old = x_new;
        x_new = iterativeFunc(x_old);

        std::cout << x_old << " " << x_new - x_old << std::endl;

        counter++;
        continueCriteria = (counter < maxIterationCount) && (std::abs(x_new - x_old) > tolerance);
    }

    return {x_new, counter < maxIterationCount};
}

#endif //COMPUTATIONALMATH_NONLINEARMPI_H
