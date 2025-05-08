#include "task7.h"

#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <vector>

#include "function.h"
#include "model/linear_ode_solver.h"
#include "util/table.h"

namespace semester6_task7 {
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
}  // namespace semester6_task7

namespace tasks {
void Semester6Task7() {
    std::cout << "Boundary value problem: proection methods" << '\n';
    std::cout << "Equation: u'' + 6x * u' - 12x^2 * u = 30x^3 - 12x^5, u(1) = 0, u(2) = 6" << '\n';
    std::cout << "Precise solution: u(x) = x^3 - x" << '\n';
    model::Func q = [](double x) { return 6 * x; };
    model::Func r = [](double x) { return -12 * x * x; };
    model::Func f = [](double x) { return 30 * x * x * x - 12 * x * x * x * x * x; };
    double a = 1;
    double b = 2;
    double a_value = 0;
    double b_value = 6;
    model::Func precise_solution = [](double x) { return x * x * x - x; };
    model::BoundaryCondition<model::BoundaryConditionKind::FirstKind> cond = {a, b, a_value,
                                                                              b_value};
    model::LinearODESolver solver({q, r}, f, cond);
    std::cout << "N = 5:" << '\n';
    auto approximate_solution = solver.Solve(5, model::ODESolveMethod::Collocation);
    semester6_task7::CompareSolutions(approximate_solution, precise_solution, a, b);

    std::cout << "Equation: u'' + 6x * u' - 12x^2 * u = 30x^3 - 12x^5, u(1) - u'(1) = -2, u(2) + "
                 "u'(2) = 17"
              << '\n';
    std::cout << "Precise solution: u(x) = x^3 - x" << '\n';
    a_value = -2;
    b_value = 17;
    model::BoundaryCondition<model::BoundaryConditionKind::ThirdKind> cond2 = {a, 1, -1, a_value,
                                                                               b, 1, 1,  b_value};
    model::LinearODESolver solver2({q, r}, f, cond2);
    std::cout << "N = 5:" << '\n';
    auto approximate_solution2 = solver2.Solve(5, model::ODESolveMethod::Collocation);
    semester6_task7::CompareSolutions(approximate_solution2, precise_solution, a, b);

    std::cout << "Equation: u'' + 6x * u' - 12x^2 * u = 30x^3 - 12x^5, u(0) = 3, u(1) = 3" << '\n';
    std::cout << "Precise solution: u(x) = 3x^6 + x^5 - 4x^4 + 3" << '\n';
    model::Func q1 = [](double x) { return -3 * x - 1; };
    model::Func r1 = [](double x) { return -x * x - 1; };
    model::Func f1 = [](double x) {
        return -3 * std::pow(x, 8) - std::pow(x, 7) - 53 * std::pow(x, 6) - 34 * std::pow(x, 5) +
               137 * std::pow(x, 4) + 36 * std::pow(x, 3) - 51 * std::pow(x, 2) - 3;
    };
    double a1 = 0;
    double b1 = 1;
    double a_value1 = 3;
    double b_value1 = 3;
    model::Func precise_solution1 = [](double x) {
        return 3 * std::pow(x, 6) + std::pow(x, 5) - 4 * std::pow(x, 4) + 3;
    };
    model::BoundaryCondition<model::BoundaryConditionKind::FirstKind> cond3 = {a1, b1, a_value1,
                                                                               b_value1};
    model::LinearODESolver solver3({q1, r1}, f1, cond3);
    std::cout << "N = 3:" << '\n';
    auto approximate_solution3 = solver3.Solve(3, model::ODESolveMethod::Collocation);
    semester6_task7::CompareSolutions(approximate_solution3, precise_solution1, a1, b1);

    std::cout << "N = 4:" << '\n';
    auto approximate_solution4 = solver3.Solve(4, model::ODESolveMethod::Collocation);
    semester6_task7::CompareSolutions(approximate_solution4, precise_solution1, a1, b1);

    std::cout << "N = 5:" << '\n';
    auto approximate_solution5 = solver3.Solve(5, model::ODESolveMethod::Collocation);
    semester6_task7::CompareSolutions(approximate_solution5, precise_solution1, a1, b1);
}

}  // namespace tasks
