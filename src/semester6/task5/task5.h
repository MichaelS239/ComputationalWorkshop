#pragma once

#include <vector>

#include "model/matrix.h"

namespace semester6_task5 {
void PrintEigen(model::Matrix const& matrix, std::vector<double> const& eps);
}  // namespace semester6_task5

namespace tasks {
void Semester6Task5();
}  // namespace tasks
