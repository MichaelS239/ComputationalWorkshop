#pragma once

#include <cstddef>
#include <functional>

#include "function.h"

namespace semester6_task7 {
void CompareSolutions(model::Func approximate_solution, model::Func galerkin_approximate_solution,
                      model::Func precise_solution, double a, double b);
}  // namespace semester6_task7

namespace tasks {
void Semester6Task7();
}  // namespace tasks
