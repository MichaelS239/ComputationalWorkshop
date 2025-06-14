#include "heat_equation_solver.h"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

#include "matrix.h"

namespace model {

std::pair<std::vector<double>, std::vector<double>> HeatEquationSolver::CreatePoints(
        std::size_t N, std::size_t K, double h, double tau) const {
    std::vector<double> x_points(N + 1);
    for (std::size_t i = 0; i != N + 1; ++i) {
        x_points[i] = i * h;
    }
    std::vector<double> t_points(K + 1);
    for (std::size_t i = 0; i != K + 1; ++i) {
        t_points[i] = i * tau;
    }

    return {std::move(x_points), std::move(t_points)};
}

TwoVariableFunc HeatEquationSolver::Solve(std::size_t N, std::size_t K,
                                          model::HeatEquationScheme scheme) const {
    double h = a / N;
    double tau = T / K;
    auto [x_points, t_points] = CreatePoints(N, K, h, tau);
    std::vector<std::vector<double>> coefs(K + 1, std::vector<double>(N + 1));
    for (std::size_t i = 0; i != N + 1; ++i) {
        coefs[0][i] = f0(x_points[i]);
    }
    for (std::size_t i = 1; i != K + 1; ++i) {
        coefs[i][0] = f1(t_points[i]);
        coefs[i][N] = f2(t_points[i]);
    }

    double gamma = tau * k / (h * h);
    if (scheme == model::HeatEquationScheme::Explicit) {
        for (std::size_t i = 1; i != K + 1; ++i) {
            for (std::size_t j = 1; j != N; ++j) {
                coefs[i][j] = gamma * coefs[i - 1][j - 1] + (1 - 2 * gamma) * coefs[i - 1][j] +
                              gamma * coefs[i - 1][j + 1] + tau * f(x_points[j], t_points[i - 1]);
            }
        }
    } else if (scheme == model::HeatEquationScheme::Implicit) {
        for (std::size_t i = 1; i != K + 1; ++i) {
            Matrix matrix(N - 1);
            std::vector<double> vec(N - 1);
            for (std::size_t j = 1; j != N; ++j) {
                if (j == 1) {
                    matrix[j - 1][j - 1] = -(1 + 2 * gamma);
                    matrix[j - 1][j] = gamma;
                    vec[j - 1] = -coefs[i - 1][j] - tau * f(x_points[j], t_points[i]) -
                                 gamma * coefs[i][j - 1];
                } else if (j == N - 1) {
                    matrix[j - 1][j - 2] = gamma;
                    matrix[j - 1][j - 1] = -(1 + 2 * gamma);
                    vec[j - 1] = -coefs[i - 1][j] - tau * f(x_points[j], t_points[i]) -
                                 gamma * coefs[i][j + 1];
                } else {
                    matrix[j - 1][j - 2] = gamma;
                    matrix[j - 1][j - 1] = -(1 + 2 * gamma);
                    matrix[j - 1][j] = gamma;
                    vec[j - 1] = -coefs[i - 1][j] - tau * f(x_points[j], t_points[i]);
                }
            }
            std::vector<double> sol = matrix.SolveTridiagonalSystem(vec);
            for (std::size_t j = 0; j != N - 1; ++j) {
                coefs[i][j + 1] = sol[j];
            }
        }
    }

    TwoVariableFunc solution = [this, x_points, t_points, coefs, h, tau](double x, double t) {
        if (x < 0 || x > a || t < 0 || t > T) return 0.0;

        std::size_t x_index =
                std::lower_bound(x_points.begin(), x_points.end(), x) - x_points.begin();
        std::size_t t_index =
                std::lower_bound(t_points.begin(), t_points.end(), t) - t_points.begin();
        return coefs[t_index][x_index];
    };

    return solution;
}

}  // namespace model
