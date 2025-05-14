#pragma once

#include <cmath>
#include <cstddef>

#include "optimization_methods.h"
#include "util/table.h"

namespace semester6_task10 {
template <std::size_t N>
void CompareSolutions(std::array<double, N> const& precise_x, double precise_value,
                      MethodInfo<N> const& approximate_solution) {
    std::vector<std::vector<double>> table(N, std::vector<double>(3));
    for (std::size_t i = 0; i != N; ++i) {
        table[i][0] = precise_x[i];
        table[i][1] = approximate_solution.local_minimum[i];
        table[i][2] = std::abs(table[i][1] - table[i][0]);
    }
    util::PrintTable(table, {"Precise solution", "Approximate solution", "Difference"});
}

}  // namespace semester6_task10

namespace tasks {
void Semester6Task10();
}  // namespace tasks
