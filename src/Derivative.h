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
struct derivativeCoeff {
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

template<typename T, indexType N>
[[nodiscard]]
Eigen::Matrix<T, N, N> fillMatrix(const std::array<T, N> &points) noexcept {

    Eigen::Matrix<T, N, N> result;

    result.col(0) = static_cast<Eigen::Matrix<scalar, N, 1>>(points.data());

    for (std::size_t i = 1; i < N; i++) {
        result.col(i) = (result.col(i - 1).asDiagonal()) * result.col(0);
    }

    return result.transpose();
}

template<typename T, indexType N, indexType derivativeOrder>
[[nodiscard]]
derivativeCoeff<T, N> calcDerivativeCoef(const std::array<T, N> &points) noexcept {

    Eigen::Vector<T, N> coeff;
    coeff(derivativeOrder - 1) = factorial(derivativeOrder);

    const Eigen::Vector<T, N> otherPointsCoeff = fillMatrix(points).colPivHouseholderQr().solve(coeff);

    return {-otherPointsCoeff.sum(), convertVectorToArray(otherPointsCoeff)};
}

#endif //COMPUTATIONALMATH_DERIVATIVE_H
