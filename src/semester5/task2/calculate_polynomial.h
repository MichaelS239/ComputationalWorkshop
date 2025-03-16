#pragma once

#include <vector>

namespace semester5_task2 {
double f(double x);
std::vector<std::pair<double, double>> CalculateTable(double a, double b, int m);
double CalculatePolynomial(std::vector<std::pair<double, double>> const& table, double x, int n);

}  // namespace semester5_task2
