#pragma once

#include <cstddef>
#include <vector>

#include "model/matrix.h"

namespace semester6_task3 {
std::pair<std::vector<std::vector<double>>, std::pair<std::size_t, std::size_t>> SolveSystem(
        model::Matrix const& matrix, double eps, std::vector<double> vector = {});
}  // namespace semester6_task3
