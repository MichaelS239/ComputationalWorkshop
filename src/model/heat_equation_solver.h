#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "function.h"
#include "solve_methods.h"

namespace model {
/**
 * Solves heat equation in the following form:
 * u_t(x, t) = k * u_xx(x, t) + f(x, t)
 * k = const > 0, 0 < x < a, 0 < t <= T
 *
 * Boundary conditions:
 * u(x, 0) = f0(x), 0 <= x <= a
 * u(0, t) = f1(t), 0 <= t <= T
 * u(a, t) = f2(t), 0 <= t <= T
 */
class HeatEquationSolver {
private:
    double k;
    double a;
    double T;
    TwoVariableFunc f;
    Func f0;
    Func f1;
    Func f2;

    std::pair<std::vector<double>, std::vector<double>> CreatePoints(std::size_t n, std::size_t k,
                                                                     double h, double tau) const;

public:
    HeatEquationSolver(double k_coef, double a_coef, double t_coef, TwoVariableFunc f_func,
                       Func f0_func, Func f1_func, Func f2_func)
        : f(f_func), f0(f0_func), f1(f1_func), f2(f2_func) {
        if (k_coef <= 0) {
            throw std::invalid_argument("Error: k must be positive");
        }
        k = k_coef;
        if (a_coef <= 0) {
            throw std::invalid_argument("Error: a must be positive");
        }
        a = a_coef;
        if (t_coef <= 0) {
            throw std::invalid_argument("Error: T must be positive");
        }
        T = t_coef;
    }

    TwoVariableFunc Solve(
            std::size_t N, std::size_t K,
            model::HeatEquationScheme scheme = model::HeatEquationScheme::Implicit) const;
};

}  // namespace model
