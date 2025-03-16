#pragma once

#include <vector>

#include "model/function.h"

namespace semester5_task3 {

double f1(double x);
double f1_prime(double x);
double f1_prime2(double x);
double f2(double x);
double f2_prime(double x);
double f2_prime2(double x);

std::vector<std::pair<double, double>> CalculateTable(double x0, double h, int m, model::Func func);
std::vector<std::vector<double>> CalculateDerivatives(
        std::vector<std::pair<double, double>> const& table, double h, model::Func func_prime,
        model::Func func_second_prime);
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> CalculateRunge(
        std::vector<std::pair<double, double>> const& table, double x0, double h, int m,
        model::Func func, model::Func func_prime, model::Func func_second_prime);

}  // namespace semester5_task3
