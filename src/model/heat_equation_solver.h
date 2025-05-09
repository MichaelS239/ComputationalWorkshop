#pragma once

#include "function.h"

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

public:
    HeatEquationSolver(double k, double a, double t, TwoVariableFunc f, Func f0, Func f1, Func f2)
        : k(k), a(a), T(t), f(f), f0(f0), f1(f1), f2(f2) {}

    TwoVariableFunc Solve(double h, double tau) const;
};

}  // namespace model
