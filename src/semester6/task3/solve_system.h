#pragma once

#include <vector>

#include "model/matrix.h"

namespace semester6_task3 {
std::vector<std::vector<double>> SolveSystem(model::Matrix const& matrix, double eps,
                                             std::vector<double> vector = {});
}  // namespace semester6_task3
