#pragma once

#include <cmath>
#include <cstddef>

#include "optimization_methods.h"
#include "util/table.h"

namespace semester6_task10 {
template <std::size_t N>
void CompareSolutions(std::array<double, N> const& precise_x, double precise_value,
                      MethodInfo<N> const& gradient_solution,
                      MethodInfo<N> const& heavy_ball_solution) {
    std::vector<std::vector<double>> table(N, std::vector<double>(5));
    for (std::size_t i = 0; i != N; ++i) {
        table[i][0] = precise_x[i];
        table[i][1] = gradient_solution.local_minimum[i];
        table[i][2] = std::abs(table[i][1] - table[i][0]);
        table[i][3] = heavy_ball_solution.local_minimum[i];
        table[i][4] = std::abs(table[i][3] - table[i][0]);
    }
    util::PrintTable(table, {"Precise solution", "Gradient descent", "Difference", "Heavy ball",
                             "Difference"});
}

}  // namespace semester6_task10

namespace tasks {
void Semester6Task10();
}  // namespace tasks
