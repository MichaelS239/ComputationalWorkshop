#pragma once

#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

#include "optimization_methods.h"
#include "util/table.h"

namespace semester6_task11 {
template <std::size_t N>
void CompareProjectionAndPenaltySolutions(std::array<double, N> const& precise_x,
                                          double precise_value, model::MultivariableFunc<N> f,
                                          model::GradientFunc<N> gradient,
                                          std::vector<Boundary<N>> const& boundaries,
                                          std::array<double, N> const& projection_x_start,
                                          std::array<double, N> const& penalty_x_start,
                                          double eps) {
    std::cout << "Gradient projection method: x0 = (";
    for (std::size_t i = 0; i != N; ++i) {
        std::cout << projection_x_start[i];
        if (i != N - 1) {
            std::cout << ", ";
        }
    }
    std::cout << ")" << '\n';
    MethodInfo gradient_projection_info =
            GradientProjectionMethod(f, gradient, boundaries, projection_x_start, eps);
    std::cout << "Number of iterations: " << gradient_projection_info.iteration_number << '\n';

    std::cout << "Penalty function method: x0 = (";
    for (std::size_t i = 0; i != N; ++i) {
        std::cout << penalty_x_start[i];
        if (i != N - 1) {
            std::cout << ", ";
        }
    }
    std::cout << ")" << '\n';
    MethodInfo penalty_function_info =
            PenaltyFunctionMethod(f, gradient, boundaries, penalty_x_start, eps);
    std::cout << "Number of iterations: " << penalty_function_info.iteration_number << '\n';

    std::vector<std::vector<double>> table(N, std::vector<double>(5));
    for (std::size_t j = 0; j != N; ++j) {
        table[j][0] = precise_x[j];
        table[j][1] = gradient_projection_info.local_minimum[j];
        table[j][2] = std::abs(table[j][1] - table[j][0]);
        table[j][3] = penalty_function_info.local_minimum[j];
        table[j][4] = std::abs(table[j][3] - table[j][0]);
    }
    util::PrintTable(table, {"Precise solution", "Gradient projection method", "Difference",
                             "Penalty function method", "Difference"});
}

template <std::size_t N>
void ComparePenaltyAndBarrierSolutions(std::array<double, N> const& precise_x, double precise_value,
                                       model::MultivariableFunc<N> f,
                                       model::GradientFunc<N> gradient,
                                       std::vector<Boundary<N>> const& boundaries,
                                       std::array<double, N> const& penalty_x_start,
                                       std::array<double, N> const& barrier_x_start, double eps) {
    std::cout << "Penalty function method: x0 = (";
    for (std::size_t i = 0; i != N; ++i) {
        std::cout << penalty_x_start[i];
        if (i != N - 1) {
            std::cout << ", ";
        }
    }
    std::cout << ")" << '\n';
    MethodInfo penalty_function_info =
            PenaltyFunctionMethod(f, gradient, boundaries, penalty_x_start, eps);
    std::cout << "Number of iterations: " << penalty_function_info.iteration_number << '\n';

    std::cout << "Barrier function method: x0 = (";
    for (std::size_t i = 0; i != N; ++i) {
        std::cout << barrier_x_start[i];
        if (i != N - 1) {
            std::cout << ", ";
        }
    }
    std::cout << ")" << '\n';
    MethodInfo barrier_function_info =
            BarrierFunctionMethod(f, gradient, boundaries, barrier_x_start, eps);
    std::cout << "Number of iterations: " << barrier_function_info.iteration_number << '\n';

    std::vector<std::vector<double>> table(N, std::vector<double>(5));
    for (std::size_t j = 0; j != N; ++j) {
        table[j][0] = precise_x[j];
        table[j][1] = penalty_function_info.local_minimum[j];
        table[j][2] = std::abs(table[j][1] - table[j][0]);
        table[j][3] = barrier_function_info.local_minimum[j];
        table[j][4] = std::abs(table[j][3] - table[j][0]);
    }
    util::PrintTable(table, {"Precise solution", "Penalty function method", "Difference",
                             "Barrier function method", "Difference"});
}

}  // namespace semester6_task11

namespace tasks {
void Semester6Task11();
}  // namespace tasks
