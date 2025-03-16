#pragma once

#include <functional>
#include <vector>

namespace semester5_task5 {
using CalcMethod = std::function<std::vector<double>(int, int, double)>;

double Solution(double x);
std::vector<double> CalculateTaylor(int start, int end, double h);
std::vector<double> CalculateAdams(int start, int end, double h);
std::vector<double> CalculateRungeKutta(int start, int end, double h);
std::vector<double> CalculateEuler(int start, int end, double h);
std::vector<double> CalculateEuler1(int start, int end, double h);
std::vector<double> CalculateEuler2(int start, int end, double h);

std::vector<std::vector<double>> CompareWithSolution(int start, int end, double h,
                                                     CalcMethod method);

}  // namespace semester5_task5
