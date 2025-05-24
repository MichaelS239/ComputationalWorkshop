#include "task8.h"

#include <cmath>
#include <iostream>

#include "model/function.h"
#include "model/heat_equation_solver.h"
#include "util/table.h"

namespace semester6_task8 {
void CompareSolutions(model::TwoVariableFunc explicit_approximate_solution,
                      model::TwoVariableFunc implicit_approximate_solution,
                      model::TwoVariableFunc precise_solution, double a, double T) {
    std::size_t points_num = 5;
    double h = a / points_num;
    double tau = T / points_num;
    std::vector<std::vector<double>> table((points_num + 1) * (points_num + 1),
                                           std::vector<double>(7));
    std::size_t index = 0;
    for (std::size_t i = 0; i != points_num + 1; ++i) {
        for (std::size_t j = 0; j != points_num + 1; ++j) {
            table[index][0] = j * h;
            table[index][1] = i * tau;
            table[index][2] = precise_solution(table[index][0], table[index][1]);
            table[index][3] = explicit_approximate_solution(table[index][0], table[index][1]);
            table[index][4] = std::abs(table[index][3] - table[index][2]);
            table[index][5] = implicit_approximate_solution(table[index][0], table[index][1]);
            table[index][6] = std::abs(table[index][5] - table[index][2]);
            ++index;
        }
    }
    util::PrintTable(table, {"x", "t", "Precise solution", "Explicit scheme", "Difference",
                             "Implicit scheme", "Difference"});
}

}  // namespace semester6_task8

namespace tasks {
void Semester6Task8() {
    std::cout << "Heat equation: finite difference methods" << '\n';
    std::cout << '\n';
    std::cout << "Equation: u_t(x, t) = 1/3 * u_xx(x, t) + 4xt - x^3, 0 <= x <= 1, 0 <= t <= 1"
              << '\n';
    std::cout << "Boundary conditions: u(x, 0) = 0, u(0, t) = 0, u(1, t) = t^2 - t" << '\n';
    std::cout << "Precise solution: u(x, t) = x * t^2 - x^3 * t" << '\n';

    double k = static_cast<double>(1) / 3;
    double a = 1;
    double T = 1;
    model::TwoVariableFunc f = [](double x, double t) { return 4 * x * t - x * x * x; };
    model::Func f0 = [](double x) { return 0.0; };
    model::Func f1 = [](double t) { return 0.0; };
    model::Func f2 = [](double t) { return t * t - t; };
    model::TwoVariableFunc precise_solution = [](double x, double t) {
        return x * t * t - x * x * x * t;
    };

    model::HeatEquationSolver solver(k, a, T, f, f0, f1, f2);
    std::cout << "h = 0.1, tau = 0.01:" << '\n';
    auto explicit_approximate_solution = solver.Solve(10, 100, model::HeatEquationScheme::Explicit);
    auto implicit_approximate_solution = solver.Solve(10, 100, model::HeatEquationScheme::Implicit);
    semester6_task8::CompareSolutions(explicit_approximate_solution, implicit_approximate_solution,
                                      precise_solution, a, T);
    std::cout << "h = 0.1, tau = 0.015:" << '\n';
    auto explicit_approximate_solution1 = solver.Solve(10, 67, model::HeatEquationScheme::Explicit);
    auto implicit_approximate_solution1 = solver.Solve(10, 67, model::HeatEquationScheme::Implicit);
    semester6_task8::CompareSolutions(explicit_approximate_solution1,
                                      implicit_approximate_solution1, precise_solution, a, T);
    std::cout << "h = 0.1, tau = 0.02:" << '\n';
    auto explicit_approximate_solution2 = solver.Solve(10, 50, model::HeatEquationScheme::Explicit);
    auto implicit_approximate_solution2 = solver.Solve(10, 50, model::HeatEquationScheme::Implicit);
    semester6_task8::CompareSolutions(explicit_approximate_solution2,
                                      implicit_approximate_solution2, precise_solution, a, T);
    std::cout << "h = 0.1, tau = 0.1:" << '\n';
    auto explicit_approximate_solution3 = solver.Solve(10, 10, model::HeatEquationScheme::Explicit);
    auto implicit_approximate_solution3 = solver.Solve(10, 10, model::HeatEquationScheme::Implicit);
    semester6_task8::CompareSolutions(explicit_approximate_solution3,
                                      implicit_approximate_solution3, precise_solution, a, T);
}

}  // namespace tasks
