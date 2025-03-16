#pragma once

#include <cmath>
#include <vector>

#include "model/function.h"

namespace semester5_task4_1 {

double f(double x);
static double const integral = -std::sin(1) + 0.946083070367183;
std::pair<std::vector<std::vector<double>>, std::vector<double>> CreateSystem(
        std::vector<double> const& points);

}  // namespace semester5_task4_1
