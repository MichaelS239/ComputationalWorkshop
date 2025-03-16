#pragma once

#include <vector>

#include "model/root_calculator.h"

namespace semester5_task1 {
double f(double x);
double f_prime(double x);
std::vector<std::pair<double, double>> FindRootSegments(double a, double b, int n);
void FindRoots(model::RootCalculator root_calculator,
               std::vector<std::pair<double, double>> root_segments, double eps);

}  // namespace semester5_task1
