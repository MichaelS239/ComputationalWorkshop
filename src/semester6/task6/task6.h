#pragma once

#include "model/function.h"
#include "model/linear_ode_solver.h"

namespace semester6_task6 {
void CompareSolutions(model::Func approximate_solution, model::Func precise_solution, double a,
                      double b);
void PrintInfo(model::SolveInfo const& info);
}  // namespace semester6_task6

namespace tasks {
void Semester6Task6();
}  // namespace tasks
