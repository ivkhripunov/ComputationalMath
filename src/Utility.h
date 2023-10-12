//
// Created by ivankhripunov on 11.10.23.
//

#ifndef COMPUTATIONALMATH_UTILITY_H
#define COMPUTATIONALMATH_UTILITY_H

#include <cmath>
#include <Eigen/Dense>

using scalar = double;
using indexType = std::size_t;
using Vector2d = Eigen::Vector2d;

struct Segment {
    scalar begin;
    scalar end;
};

#endif //COMPUTATIONALMATH_UTILITY_H
