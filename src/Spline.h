//
// Created by ivankhripunov on 16.10.23.
//

#ifndef COMPUTATIONALMATH_SPLINE_H
#define COMPUTATIONALMATH_SPLINE_H

#include <vector>
#include <type_traits>
#include <iostream>
#include <initializer_list>
#include "Utility.h"

template<typename Type>
class Trio {
public:
    Type first_element;
    Type second_element;
    Type third_element;

    Trio(const Type &first_value, const Type &second_value, const Type &third_value) : first_element(first_value),
                                                                                       second_element(second_value),
                                                                                       third_element(third_value) {}

    [[nodiscard]] Type operator[](std::size_t i) const {
        switch (i) {
            case 0:
                return first_element;
            case 1:
                return second_element;
            case 2:
                return third_element;
            default:
                break;
        }

    }

    [[nodiscard]] Type &operator[](std::size_t i) {
        switch (i) {
            case 0:
                return first_element;
            case 1:
                return second_element;
            case 2:
                return third_element;
            default:
                break;
        }

    }

};

template<typename Type>
class TridiagonalMatrix {
private:
    std::vector<Trio<Type>> matrix;
    indexType size_;

public:
    TridiagonalMatrix(const std::vector<Type> &lower_diagonal,
                      const std::vector<Type> &main_diagonal,
                      const std::vector<Type> &upper_diagonal) {

        size_ = main_diagonal.size();

        matrix.reserve(size_);

        for (std::size_t i = 0; i < size_; ++i) {
            Trio tmp(lower_diagonal[i], main_diagonal[i], upper_diagonal[i]);
            matrix.push_back(tmp);
        }
    }

    TridiagonalMatrix(const std::initializer_list<Type> &init_lower_diagonal,
                      const std::initializer_list<Type> &init_main_diagonal,
                      const std::initializer_list<Type> &init_upper_diagonal) {

        size_ = init_main_diagonal.size();

        matrix.reserve(size_);

        for (std::size_t i = 0; i < size_; ++i) {
            Trio tmp(*(init_lower_diagonal.begin() + i),
                     *(init_main_diagonal.begin() + i),
                     *(init_upper_diagonal.begin() + i));
            matrix[i] = tmp;
        }

    }

    [[nodiscard]] size_t size() const {
        return size_;
    }

    [[nodiscard]] Trio<Type> &operator[](std::size_t i) {
        return matrix[i];
    }

    [[nodiscard]] Trio<Type> operator[](std::size_t i) const {
        return matrix[i];
    }

    friend std::ostream &operator<<(std::ostream &os, const TridiagonalMatrix<Type> &Matrix) {
        for (size_t i = 0; i < Matrix.size_; ++i) {
            std::cout << Matrix.matrix[i].first_element << " "
                      << Matrix.matrix[i].second_element << " "
                      << Matrix.matrix[i].third_element << std::endl;
        }

        return os;
    }
};


template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

template<typename mType, typename cType>
void calc_coeffs(const TridiagonalMatrix<mType> &Matrix, const std::vector<cType> &right_hand_column,
                 std::vector<DivisType<cType, mType>> &p_coeffs, std::vector<DivisType<cType, mType>> &q_coeffs) {

    p_coeffs[0] = -Matrix[0].third_element / Matrix[0].second_element;
    q_coeffs[0] = right_hand_column[0] / Matrix[0].second_element;

    for (size_t i = 0; i < Matrix.size() - 1; ++i) {
        p_coeffs[i + 1] =
                -Matrix[i + 1].third_element /
                (Matrix[i + 1].first_element * p_coeffs[i] + Matrix[i + 1].second_element);

        q_coeffs[i + 1] =
                (right_hand_column[i + 1] - Matrix[i + 1].first_element * q_coeffs[i]) /
                (Matrix[i + 1].first_element * p_coeffs[i] + Matrix[i + 1].second_element);

    }
}

/** Функция для решения методм  прогонки **/
template<typename mType, typename cType>
std::vector<DivisType<cType, mType>> solve(const TridiagonalMatrix<mType> &matrix,
                                           const std::vector<cType> &column) {

    size_t size_ = matrix.size();
    std::vector<DivisType<cType, mType>> result(size_);

    std::vector<DivisType<cType, mType>> p_values(size_);

    calc_coeffs(matrix, column, p_values, result);

    for (int i = size_ - 2; i >= 0; --i) {
        result[i] = p_values[i] * result[i + 1] + result[i];
    }
    return result;
}

template<typename T>
struct SplinePolynom {
    T a;  //zero degree
    T b;  //first degree
    T c;  //second degree
    T d;  //third degree
};


template<typename xType, typename yType>
class CubicSpline {

    std::vector<xType> points_;
    std::vector<SplinePolynom<DivisType<xType, yType>>> splines_;

public:
    CubicSpline(const std::vector<xType> &points,
                const std::vector<yType> &values) : points_(points) {

        const indexType pointsCount = points.size();
        const indexType splinesCount = pointsCount - 1;
        const indexType eqCount = pointsCount - 2;
        const indexType sideDiagonalSize = pointsCount - 3;

        splines_.resize(splinesCount);

        std::vector<DivisType<xType, yType>>
                h_array(splinesCount),
                uDouble_array(splinesCount),
                rightHandSide(eqCount);

        for (indexType i = 0; i < splinesCount; ++i) {
            h_array[i] = points[i + 1] - points[i];
            uDouble_array[i] = (values[i + 1] - values[i]) / h_array[i];
        }

        for (indexType i = 0; i < eqCount; ++i) {
            rightHandSide[i] = 6 * (uDouble_array[i + 1] - uDouble_array[i]) / (h_array[i + 1] + h_array[i]);
        }

        const std::vector<DivisType<xType, yType>> mainDiagonal(eqCount, 2);
        std::vector<DivisType<xType, yType>> lowerDiagonal(eqCount), upperDiagonal(eqCount);

        lowerDiagonal[0] = 0;
        upperDiagonal[eqCount - 1] = 0;
        for (indexType i = 0; i < sideDiagonalSize; ++i) {
            const DivisType<xType, yType> denominator = h_array[i + 1] + h_array[i];
            upperDiagonal[i] = h_array[i + 1] / denominator;
            lowerDiagonal[i + 1] = h_array[i] / denominator;
        }

        const TridiagonalMatrix<xType> matrix(lowerDiagonal, mainDiagonal, upperDiagonal);

        // c_0 ... c_N-2
        const std::vector<DivisType<xType, yType>> cResult = solve(matrix, rightHandSide);


        //set c
        splines_[splinesCount - 1].c = 0;
        for (indexType i = 0; i < splinesCount - 1; ++i) {
            splines_[i].c = cResult[i];
        }


        //set d
        splines_[0].d = splines_[0].c / h_array[0];
        for (indexType i = 1; i < splinesCount; ++i) {
            splines_[i].d = (splines_[i].c - splines_[i - 1].c) / h_array[i];
        }


        //set a
        for (indexType i = 0; i < splinesCount; ++i) {
            splines_[i].a = values[i + 1];
        }


        //set b
        splines_[0].b = splines_[0].c * h_array[0] / 3 + uDouble_array[0];
        for (indexType i = 1; i < splinesCount; ++i) {
            splines_[i].b = (2 * splines_[i].c + splines_[i - 1].c) * h_array[i] / 6 + uDouble_array[i];
        }


//        for (indexType j = 1; j < splinesCount; ++j) {
//            std::cout << splines_[j].a
//                         - splines_[j].b * h_array[j]
//                         + splines_[j].c * h_array[j] * h_array[j] / 2
//                         - splines_[j].d * h_array[j] * h_array[j] * h_array[j] / 6
//                         - splines_[j - 1].a << std::endl;
//
//            std::cout << splines_[j - 1].b
//                         + splines_[j].c * h_array[j]
//                         - splines_[j].d * h_array[j] * h_array[j] / 2
//                         - splines_[j].b << std::endl;
//
//        }

    };

    [[nodiscard]]
    yType interpolate(const xType &x) const noexcept {

        for (indexType i = 0; i < splines_.size(); ++i) {

            if (points_[i] <= x and x <= points_[i + 1]) {
                const xType delta = x - points_[i + 1];

                return ((splines_[i].d / 6 * delta + splines_[i].c / 2) * delta + splines_[i].b) *
                       delta + splines_[i].a;

//                return splines_[i].a + splines_[i].b * delta
//                       + splines_[i].c * delta * delta / 2
//                       + splines_[i].d * delta * delta * delta / 6;

            }
        }
    };

    [[nodiscard]]
    yType extrapolateRight(const xType &x) const noexcept {

        const indexType splinesCount = splines_.size();
        const indexType rightSplineIndex = splinesCount - 1;
        const xType delta = x - points_[points_.size() - 1];

        return ((splines_[rightSplineIndex].d / 6 * delta + splines_[rightSplineIndex].c / 2) * delta + splines_[rightSplineIndex].b) *
               delta + splines_[rightSplineIndex].a;
    }
};

#endif //COMPUTATIONALMATH_SPLINE_H
