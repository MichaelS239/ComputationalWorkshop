#pragma once

#include "model/function.h"

namespace semester6_task8 {
void CompareSolutions(model::TwoVariableFunc explicit_approximate_solution,
                      model::TwoVariableFunc implicit_approximate_solution,
                      model::TwoVariableFunc precise_solution, double a, double T);
}  // namespace semester6_task8

namespace tasks {
void Semester6Task8();
}  // namespace tasks
