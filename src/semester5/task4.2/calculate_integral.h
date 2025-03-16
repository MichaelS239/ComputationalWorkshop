#pragma once

#include <string>
#include <vector>

#include "model/function.h"

namespace semester5_task4_2 {

double f(double x);
double Integral(double a, double b);
std::vector<std::pair<std::string, double>> CalculateIntegrals(double a, double b, model::Func f);

}  // namespace semester5_task4_2
