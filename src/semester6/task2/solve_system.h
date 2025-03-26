#pragma once

#include <vector>

#include "model/matrix.h"

namespace semester6_task2 {
std::vector<std::vector<double>> SolveSystem(model::Matrix const& matrix,
                                             std::vector<double> vector = {});
}  // namespace semester6_task2
