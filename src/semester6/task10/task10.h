#pragma once

#include <cmath>
#include <cstddef>
#include <iostream>

#include "optimization_methods.h"
#include "util/table.h"

namespace semester6_task10 {
template <std::size_t N>
void CompareSolutions(std::array<double, N> const& precise_x, double precise_value,
                      model::MultivariableFunc<N> f, model::GradientFunc<N> gradient,
                      model::SecondDerivativeFunc<N> second_derivative, double eps,
                      std::vector<std::pair<double, double>> const& coefs) {
    std::vector<std::vector<double>> iter_table(coefs.size(), std::vector<double>(6));

    for (std::size_t i = 0; i != coefs.size(); ++i) {
        double alpha = coefs[i].first;
        double beta = coefs[i].second;
        MethodInfo gradient_info = GradientDescent(f, gradient, eps, alpha);
        MethodInfo heavy_ball_info = HeavyBallMethod(f, gradient, eps, alpha, beta);
        MethodInfo nesterov_info = NesterovMethod(f, gradient, eps, alpha);
        MethodInfo newton_info = NewtonMethod(f, gradient, second_derivative, eps);

        iter_table[i] = {alpha,
                         beta,
                         static_cast<double>(gradient_info.iteration_number),
                         static_cast<double>(heavy_ball_info.iteration_number),
                         static_cast<double>(nesterov_info.iteration_number),
                         static_cast<double>(newton_info.iteration_number)};

        std::vector<std::vector<double>> table(N, std::vector<double>(9));
        for (std::size_t j = 0; j != N; ++j) {
            table[j][0] = precise_x[j];
            table[j][1] = gradient_info.local_minimum[j];
            table[j][2] = std::abs(table[j][1] - table[j][0]);
            table[j][3] = heavy_ball_info.local_minimum[j];
            table[j][4] = std::abs(table[j][3] - table[j][0]);
            table[j][5] = nesterov_info.local_minimum[j];
            table[j][6] = std::abs(table[j][5] - table[j][0]);
            table[j][7] = newton_info.local_minimum[j];
            table[j][8] = std::abs(table[j][7] - table[j][0]);
        }
        std::cout << "Step: " << alpha << ", momentum: " << beta << '\n';
        util::PrintTable(table, {"Precise solution", "Gradient descent", "Difference",
                                 "Heavy ball method", "Difference", "Nesterov method", "Difference",
                                 "Newton method", "Difference"});
    }
    std::cout << "Iteration numbers:" << '\n';
    util::PrintTable(iter_table, {"Step", "Momentum", "Gradient descent", "Heavy ball method",
                                  "Nesterov method", "Newton method"});
}

}  // namespace semester6_task10

namespace tasks {
void Semester6Task10();
}  // namespace tasks
