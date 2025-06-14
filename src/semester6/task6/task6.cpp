#include "task6.h"

#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <vector>

#include "model/function.h"
#include "model/linear_ode_solver.h"
#include "util/table.h"

namespace semester6_task6 {
void CompareSolutions(model::Func approximate_solution, model::Func precise_solution, double a,
                      double b) {
    std::size_t points_num = 10;
    double h = (b - a) / points_num;
    std::vector<std::vector<double>> table(points_num + 1, std::vector<double>(4));
    for (std::size_t i = 0; i != points_num + 1; ++i) {
        table[i][0] = a + i * h;
        table[i][1] = precise_solution(a + i * h);
        table[i][2] = approximate_solution(a + i * h);
        table[i][3] = std::abs(table[i][2] - table[i][1]);
    }
    util::PrintTable(table, {"Point", "Precise solution", "Approximate solution", "Difference"});
}

void PrintInfo(model::SolveInfo const& info) {
    std::vector<std::vector<double>> table(info.size(), std::vector<double>(4));
    for (std::size_t i = 0; i != info.size(); ++i) {
        if (i == 0) {
            table[i] = {static_cast<double>(i + 1), -log10(info[i].first), log10(info[i].second),
                        0};
        } else {
            table[i] = {static_cast<double>(i + 1), -log10(info[i].first), log10(info[i].second),
                        (log10(info[i].second) - log10(info[i - 1].second)) /
                                (-log10(info[i].first) + log10(info[i - 1].first))};
        }
    }
    util::PrintTable(table, {"Iteration number", "log(1/h)", "log(Norm)", "k"});
}

}  // namespace semester6_task6

namespace tasks {
void Semester6Task6() {
    std::cout << "Boundary value problem: finite difference method" << '\n';
    std::cout << '\n';
    std::cout << "Equation: u'' + 6x * u' - 12x^2 * u = 30x^3 - 12x^5, u(1) = 0, u(2) = 6" << '\n';
    std::cout << "Precise solution: u(x) = x^3 - x" << '\n';
    model::Func p = [](double x) { return 6 * x; };
    model::Func q = [](double x) { return -12 * x * x; };
    model::Func f = [](double x) { return 30 * x * x * x - 12 * x * x * x * x * x; };
    double a = 1;
    double b = 2;
    double a_value = 0;
    double b_value = 6;
    model::Func precise_solution = [](double x) { return x * x * x - x; };
    model::BoundaryCondition<model::BoundaryConditionKind::FirstKind> cond = {a, b, a_value,
                                                                              b_value};
    model::LinearODESolver solver({p, q}, f, cond);
    std::cout << "Precision: 1e-2" << '\n';
    auto [approximate_solution, info] = solver.Solve(1e-2);
    semester6_task6::PrintInfo(info);
    semester6_task6::CompareSolutions(approximate_solution, precise_solution, a, b);
    std::cout << "Precision: 1e-5" << '\n';
    auto [approximate_solution1, info1] = solver.Solve(1e-5);
    semester6_task6::PrintInfo(info1);
    semester6_task6::CompareSolutions(approximate_solution1, precise_solution, a, b);
    std::cout << "Precision: 1e-7" << '\n';
    auto [approximate_solution2, info2] = solver.Solve(1e-7);
    semester6_task6::PrintInfo(info2);
    semester6_task6::CompareSolutions(approximate_solution2, precise_solution, a, b);

    std::cout << "Equation: u'' + 6x * u' - 12x^2 * u = 30x^3 - 12x^5, u'(1) = 2, u'(2) = 11"
              << '\n';
    std::cout << "Precise solution: u(x) = x^3 - x" << '\n';
    a_value = 2;
    b_value = 11;
    model::BoundaryCondition<model::BoundaryConditionKind::SecondKind> cond1 = {a, b, a_value,
                                                                                b_value};
    model::LinearODESolver solver1({p, q}, f, cond1);
    std::cout << "Precision: 1e-5" << '\n';
    auto [approximate_solution3, info3] = solver1.Solve(1e-5);
    semester6_task6::PrintInfo(info3);
    semester6_task6::CompareSolutions(approximate_solution3, precise_solution, a, b);

    std::cout << "Equation: u'' + 6x * u' - 12x^2 * u = 30x^3 - 12x^5, u(1) - u'(1) = -2, u(2) + "
                 "u'(2) = 17"
              << '\n';
    std::cout << "Precise solution: u(x) = x^3 - x" << '\n';
    a_value = -2;
    b_value = 17;
    model::BoundaryCondition<model::BoundaryConditionKind::ThirdKind> cond2 = {a, 1, -1, a_value,
                                                                               b, 1, 1,  b_value};
    model::LinearODESolver solver2({p, q}, f, cond2);
    std::cout << "Precision: 1e-5" << '\n';
    auto [approximate_solution4, info4] = solver2.Solve(1e-5);
    semester6_task6::PrintInfo(info4);
    semester6_task6::CompareSolutions(approximate_solution4, precise_solution, a, b);
}

}  // namespace tasks
