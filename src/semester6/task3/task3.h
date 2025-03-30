#pragma once

#include <vector>

#include "model/matrix.h"

namespace semester6_task3 {
void PrintSystem(model::Matrix const& matrix, std::vector<double> const& eps,
                 std::vector<double> vector = {});
}  // namespace semester6_task3

namespace tasks {
void Semester6Task3();
}  // namespace tasks
