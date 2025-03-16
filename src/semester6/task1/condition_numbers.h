#pragma once

#include <cstddef>
#include <vector>

#include "model/matrix.h"

namespace semester6_task1 {
std::pair<std::vector<double>, std::vector<std::vector<double>>> CalculateConds(
        model::Matrix const& matrix);

}  // namespace semester6_task1
