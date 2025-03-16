#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

#include "model/function.h"
#include "model/polynomial.h"

namespace semester5_task6 {

double f(double x);
static double const integral = -std::sin(1) + 0.946083070367183;

std::vector<double> CalculateMoments(std::size_t n);
std::vector<double> CalculateGaussMoments(std::size_t n);
std::vector<double> CalculateMellerMoments(std::size_t n);
std::pair<std::vector<std::vector<double>>, std::vector<double>> CreatePointsSystem(
        std::size_t num_points, std::vector<double> const& moments);

}  // namespace semester5_task6
