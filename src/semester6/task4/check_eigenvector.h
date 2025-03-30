#pragma once

#include <vector>

#include "model/matrix.h"

namespace semester6_task4 {
std::vector<std::vector<double>> CheckEigenvector(model::Matrix const& matrix,
                                                  std::vector<double> const& eigenvector);
}  // namespace semester6_task4
