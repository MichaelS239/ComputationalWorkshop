#pragma once

#include <functional>

#include "model/linear_ode_solver.h"

namespace semester6_task6 {
void CompareSolutions(std::function<double(double)> approximate_solution,
                      std::function<double(double)> precise_solution, double a, double b);
void PrintInfo(model::SolveInfo const& info);
}  // namespace semester6_task6

namespace tasks {
void Semester6Task6();
}  // namespace tasks
