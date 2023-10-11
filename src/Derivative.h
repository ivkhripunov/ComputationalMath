//
// Created by ivankhripunov on 11.10.23.
//

#ifndef COMPUTATIONALMATH_DERIVATIVE_H
#define COMPUTATIONALMATH_DERIVATIVE_H

#include "Utility.h"
#include <array>

[[nodiscard]] constexpr indexType factorial(indexType n) noexcept {
    return n > 0 ? n * factorial(n - 1) : 1;
}

template<typename T, indexType N>
struct calcDerivativeCoeff {
    T centralPointCoeff;
    std::array<T, N> otherPointsCoeff;
};

template<typename T, int N>
[[nodiscard]]
std::array<T, N> convertVectorToArray(const Eigen::Vector<T, N> &vector) noexcept {
    std::array<T, N> result;

    for (indexType i = 0; i < N; ++i) {
        result[i] = vector[i];
    }
    return result;
}

template<typename T, indexType N, indexType derivativeOrder>
[[nodiscard]]
calcDerivativeCoeff<T, N> calcDerivativeCoef(const std::array<T, N> &points) noexcept {
    Eigen::Matrix<T, N, N> matrix;

    Eigen::Vector<T, N> coeff;
    coeff(derivativeOrder - 1) = factorial(derivativeOrder);

    matrix.col(0) = static_cast<Eigen::Matrix<scalar, N, 1>>(points.data());

    for (std::size_t i = 1; i < N; i++) {
        matrix.col(i) = (matrix.col(i - 1).asDiagonal()) * matrix.col(0);
    }

    const Eigen::Vector<T, N> otherPointsCoeff = matrix.transpose().colPivHouseholderQr().solve(coeff);

    return {-otherPointsCoeff.sum(), convertVectorToArray(otherPointsCoeff)};
}

#endif //COMPUTATIONALMATH_DERIVATIVE_H
