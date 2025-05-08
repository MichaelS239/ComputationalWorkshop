#pragma once

#include <cstddef>
#include <functional>

namespace semester6_task7 {
void CompareSolutions(std::function<double(double)> approximate_solution,
                      std::function<double(double)> precise_solution, double a, double b);
}  // namespace semester6_task7

namespace tasks {
void Semester6Task7();
}  // namespace tasks
