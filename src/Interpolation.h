//
// Created by ivankhripunov on 12.10.23.
//

#ifndef COMPUTATIONALMATH_INTERPOLATION_H
#define COMPUTATIONALMATH_INTERPOLATION_H

#include <array>
#include "Utility.h"

template<typename xType, typename yType, indexType N>
class NewtonInterpolator {
private:
    std::array<xType, N> xPoints_;
    std::array<yType, N> dividedDifferences_;

public:
    NewtonInterpolator(const std::array<xType, N> &xPoints, const std::array<yType, N> &yPoints) noexcept
            : xPoints_{xPoints}, dividedDifferences_(yPoints) {

        for (indexType i = 0; i < N - 1; i++) {
            for (indexType j = N - 1; j > i; j--) {
                dividedDifferences_[j] =
                        (dividedDifferences_[j] - dividedDifferences_[j - 1]) / (xPoints[j] - xPoints[j - 1 - i]);
            }
        }

        std::cout << dividedDifferences_ << std::endl << xPoints_ << std::endl;
    };

    [[nodiscard]]
    yType interpolate(const xType &x) const noexcept {

        yType result = dividedDifferences_[N - 1];
        for (std::size_t i = N - 1; i > 0; i--) {
            //суммируем горнером, добавим fma, чтобы была одинарная арифметическая ошибка
            // <<==>> result = dividedDifferences_[i - 1] + result * (x - xPoints_[i - 1]);
            result = std::fma(result, (x - xPoints_[i - 1]), dividedDifferences_[i - 1]);
        }
        return result;
    };

};

#endif //COMPUTATIONALMATH_INTERPOLATION_H
