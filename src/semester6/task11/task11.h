#pragma once

#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

#include "optimization_methods.h"
#include "util/table.h"

namespace semester6_task11 {
template <std::size_t N>
void CompareSolutions(
        std::array<double, N> const& precise_x, double precise_value, model::MultivariableFunc<N> f,
        model::GradientFunc<N> gradient, std::vector<model::MultivariableFunc<N>> const& boundaries,
        std::vector<std::vector<model::MultivariableFunc<N>>> const& boundaries_derivatives,
        std::array<double, N> x_start, double eps) {
    MethodInfo gradient_projection_info =
            GradientProjectionMethod(f, gradient, boundaries, boundaries_derivatives, x_start, eps);
    std::cout << "Number of iterations: " << gradient_projection_info.iteration_number << '\n';

    std::vector<std::vector<double>> table(N, std::vector<double>(3));
    for (std::size_t j = 0; j != N; ++j) {
        table[j][0] = precise_x[j];
        table[j][1] = gradient_projection_info.local_minimum[j];
        table[j][2] = std::abs(table[j][1] - table[j][0]);
    }
    util::PrintTable(table, {"Precise solution", "Gradient projection method", "Difference"});
}

}  // namespace semester6_task11

namespace tasks {
void Semester6Task11();
}  // namespace tasks
