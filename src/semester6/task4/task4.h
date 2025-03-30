#pragma once

#include <vector>

#include "model/matrix.h"

namespace semester6_task4 {
void PrintEigen(model::Matrix const& matrix, std::vector<double> const& eps);
}  // namespace semester6_task4

namespace tasks {
void Semester6Task4();
}  // namespace tasks
