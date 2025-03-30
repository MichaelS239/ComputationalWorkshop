#pragma once

#include <vector>

#include "model/matrix.h"

namespace semester6_task4 {
struct EigenInfo {
    std::vector<std::vector<double>> power_table;
    std::vector<std::vector<double>> scalar_table;
    double power_eigenvalue;
    double scalar_eigenvalue;
    std::size_t power_iter_num;
    std::size_t scalar_iter_num;
};

EigenInfo CalculateEigenvector(model::Matrix const& matrix, double eps);
}  // namespace semester6_task4
