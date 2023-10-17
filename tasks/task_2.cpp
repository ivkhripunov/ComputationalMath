//
// Created by ivankhripunov on 17.10.23.
//

#include <gtest/gtest.h>
#include "Integrate.h"
#include "Interpolation.h"

TEST(TASK2, TRAPSIMPGAUSS) {

    const auto function = [](const scalar x) { return std::sin(100 * x) * std::exp(-x * x) * std::cos(2 * x); };

    const Segment segment = {0, 3};

    const indexType N0 = 250;
    const std::vector<scalar> grid0 = calcGrid(N0, segment);

    const indexType N1 = 500;
    const std::vector<scalar> grid1 = calcGrid(N1, segment);

    const indexType N2 = 1000;
    const std::vector<scalar> grid2 = calcGrid(N2, segment);

    const scalar referenceResult = 0.010006097860332;

    // grid 250
    const scalar resultTrapecia250 = integrateGridTrapecia(function, grid0);
    const scalar resultSimpson250 = integrateGridSimpson(function, grid0);

    std::cout << std::endl << "#Grid 250 points:" << std::endl;

    std::cout << "Трапеция: " << resultTrapecia250 << std::endl << "Симпсон: " << resultSimpson250 << " " << std::endl;

    // grid 500
    const scalar resultTrapecia500 = integrateGridTrapecia(function, grid1);
    const scalar resultSimpson500 = integrateGridSimpson(function, grid1);

    std::cout << std::endl << "#Grid 500 points:" << std::endl;

    std::cout << "Трапеция: " << resultTrapecia500 << std::endl << "Симпсон: " << resultSimpson500 << " " << std::endl;

    //grid 1000
    const scalar resultTrapecia1000 = integrateGridTrapecia(function, grid2);
    const scalar resultSimpson1000 = integrateGridSimpson(function, grid2);

    std::cout << std::endl << "#Grid 1000 points:" << std::endl;

    std::cout << "Трапеция: " << resultTrapecia1000 << std::endl << "Симпсон: " << resultSimpson1000 << " "
              << std::endl;

    //Порядки сходимости

    const scalar p_trapecia = std::log2(
            (resultTrapecia250 - resultTrapecia500) / (resultTrapecia500 - resultSimpson1000));

    std::cout << std::endl << "Порядок сходимости трапеции (2): " << p_trapecia << std::endl;

    const scalar p_simpson = std::log2(
            (resultSimpson250 - resultSimpson500) / (resultSimpson500 - resultSimpson1000));

    std::cout << "Порядок сходимости Cимпмона (4): " << p_simpson << std::endl;

    //Экстраполяция Ричардсона

    const scalar epsTrapecia500 = (resultTrapecia250 - resultTrapecia500) / (pow(2, p_trapecia) - 1);
    const scalar epsTrapecia1000 = (resultTrapecia500 - resultTrapecia1000) / (pow(2, p_trapecia) - 1);

    const scalar resultTrapecia500_rich = resultTrapecia500 - epsTrapecia500;
    const scalar resultTrapecia1000_rich = resultTrapecia1000 - epsTrapecia1000;


    const scalar epsSimpson500 = (resultSimpson250 - resultSimpson500) / (pow(2, p_simpson) - 1);
    const scalar epsSimpson1000 = (resultSimpson500 - resultSimpson1000) / (pow(2, p_simpson) - 1);

    const scalar resultSimpson500_rich = resultSimpson500 - epsSimpson500;
    const scalar resultSimpson1000_rich = resultSimpson1000 - epsSimpson1000;

    std::cout << std::endl << "После экстраполяции Ричардсона: " << std::endl;

    std::cout << std::endl << "#Grid 500 points:" << std::endl;

    std::cout << "Трапеция: " << resultTrapecia500_rich << std::endl << "Симпсон: " << resultSimpson500_rich << " "
              << std::endl;

    std::cout << std::endl << "#Grid 1000 points:" << std::endl;

    std::cout << "Трапеция: " << resultTrapecia1000_rich << std::endl << "Симпсон: " << resultSimpson1000_rich << " "
              << std::endl;

}

TEST(TASK2, IMPROVEMENT) {

    const auto functionForInterpolation = [](const scalar x) { return std::exp(-x * x); };

    const Segment segment = {0, 3};
    const indexType N = 3;
    const std::vector<scalar> grid = calcGrid(N, segment);

    std::array<scalar, N> gridArray{};
    std::array<scalar, N> valuesArray{};
    for (indexType i = 0; i < N; ++i) {
        gridArray[i] = grid[i];
        valuesArray[i] = functionForInterpolation(grid[i]);
    }

    const NewtonInterpolator<scalar, scalar, N> newtonInterpolator(gridArray, valuesArray);

    const scalar a = 0.175406;
    const scalar b = -0.85951;
    const scalar c = 1;

    const auto integral = [&](const scalar x) {

        const scalar cos98 = std::cos(98 * x);
        const scalar cos102 = std::cos(102 * x);

        const scalar sin98 = std::sin(98 * x);
        const scalar sin102 = std::sin(102 * x);

        const scalar xSqr = x * x;

        const scalar first = -cos98 * (a * (4802 * xSqr - 1) + 4802 * (b * x + c)) / 941192;
        const scalar second = -cos102 * (a * (5202 * xSqr - 1) + 5202 * (b * x + c)) / 1061208;
        const scalar third = (2601 * sin98 + 2401 * sin102) * (2 * a * x + b) / 49960008;

        return first + second + third;
    };

    std::cout << "Значение интеграла:" << std::endl;
    std::cout << integral(3) - integral(0);
}

