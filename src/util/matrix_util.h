#pragma once

#include <vector>

namespace util {
std::vector<double> SolveSystem(std::vector<std::vector<double>> const& matrix,
                                std::vector<double> const& vector);
std::vector<std::vector<double>> CalculateInverseMatrix(
        std::vector<std::vector<double>> const& matrix);
}  // namespace util
