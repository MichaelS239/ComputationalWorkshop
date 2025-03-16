#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "model/function.h"

namespace semester5_task4_3 {
static std::vector<std::size_t> const precisions = {1, 1, 2, 2, 4};

double f(double x);
double Integral(double a, double b);
std::vector<std::pair<std::string, double>> CalculateIntegrals(double a, double b, int n,
                                                               model::Func f);
std::vector<std::pair<std::string, double>> CalculateRunge(
        std::vector<std::pair<std::string, double>> const& integrals,
        std::vector<std::pair<std::string, double>> const& new_integrals, int l);

}  // namespace semester5_task4_3
