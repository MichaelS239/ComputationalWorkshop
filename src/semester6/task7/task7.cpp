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
void CompareSolutions(model::Func collocation_approximate_solution,
                      model::Func galerkin_approximate_solution, model::Func precise_solution,
                      double a, double b) {
    std::size_t points_num = 10;
    double h = (b - a) / points_num;
    std::vector<std::vector<double>> table(points_num + 1, std::vector<double>(6));
    for (std::size_t i = 0; i != points_num + 1; ++i) {
        table[i][0] = a + i * h;
        table[i][1] = precise_solution(a + i * h);
        table[i][2] = collocation_approximate_solution(a + i * h);
        table[i][3] = std::abs(table[i][2] - table[i][1]);
        table[i][4] = galerkin_approximate_solution(a + i * h);
        table[i][5] = std::abs(table[i][4] - table[i][1]);
    }
    util::PrintTable(table, {"Point", "Precise solution", "Collocation method", "Difference",
                             "Galerkin method", "Difference"});
}
}  // namespace semester6_task7

namespace tasks {
void Semester6Task7() {
    std::cout << "Boundary value problem: proection methods" << '\n';
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
    std::cout << "N = 5:" << '\n';
    auto collocation_approximate_solution = solver.Solve(5, model::ODESolveMethod::Collocation);
    auto galerkin_approximate_solution = solver.Solve(5, model::ODESolveMethod::Galerkin);
    semester6_task7::CompareSolutions(collocation_approximate_solution,
                                      galerkin_approximate_solution, precise_solution, a, b);

    std::cout << "Equation: u'' + 6x * u' - 12x^2 * u = 30x^3 - 12x^5, u(1) - u'(1) = -2, u(2) + "
                 "u'(2) = 17"
              << '\n';
    std::cout << "Precise solution: u(x) = x^3 - x" << '\n';
    a_value = -2;
    b_value = 17;
    model::BoundaryCondition<model::BoundaryConditionKind::ThirdKind> cond2 = {a, 1, -1, a_value,
                                                                               b, 1, 1,  b_value};
    model::LinearODESolver solver2({p, q}, f, cond2);
    std::cout << "N = 5:" << '\n';
    auto collocation_approximate_solution2 = solver2.Solve(5, model::ODESolveMethod::Collocation);
    auto galerkin_approximate_solution2 = solver2.Solve(5, model::ODESolveMethod::Galerkin);
    semester6_task7::CompareSolutions(collocation_approximate_solution2,
                                      galerkin_approximate_solution2, precise_solution, a, b);

    std::cout << "Equation: u'' + 6x * u' - 12x^2 * u = 30x^3 - 12x^5, u(0) = 3, u(1) = 3" << '\n';
    std::cout << "Precise solution: u(x) = 3x^6 + x^5 - 4x^4 + 3" << '\n';
    model::Func p1 = [](double x) { return -3 * x - 1; };
    model::Func q1 = [](double x) { return -x * x - 1; };
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
    model::LinearODESolver solver3({p1, q1}, f1, cond3);
    std::cout << "N = 3:" << '\n';
    auto collocation_approximate_solution3 = solver3.Solve(3, model::ODESolveMethod::Collocation);
    auto galerkin_approximate_solution3 = solver3.Solve(3, model::ODESolveMethod::Galerkin);
    semester6_task7::CompareSolutions(collocation_approximate_solution3,
                                      galerkin_approximate_solution3, precise_solution1, a1, b1);

    std::cout << "N = 4:" << '\n';
    auto collocation_approximate_solution4 = solver3.Solve(4, model::ODESolveMethod::Collocation);
    auto galerkin_approximate_solution4 = solver3.Solve(4, model::ODESolveMethod::Galerkin);
    semester6_task7::CompareSolutions(collocation_approximate_solution4,
                                      galerkin_approximate_solution4, precise_solution1, a1, b1);

    std::cout << "N = 5:" << '\n';
    auto collocation_approximate_solution5 = solver3.Solve(10, model::ODESolveMethod::Collocation);
    auto galerkin_approximate_solution5 = solver3.Solve(10, model::ODESolveMethod::Galerkin);
    semester6_task7::CompareSolutions(collocation_approximate_solution5,
                                      galerkin_approximate_solution5, precise_solution1, a1, b1);

    std::cout << "Equation: u'' + x^2 * u' - x * u = 6 / x^4 - 3 / x, u(1) = 1, 3u(2) + u'(2) = 0.5"
              << '\n';
    std::cout << "Precise solution: u(x) = 1 / x^2" << '\n';
    model::Func p2 = [](double x) { return x * x; };
    model::Func q2 = [](double x) { return -x; };
    model::Func f2 = [](double x) { return 6 / std::pow(x, 4) - 3 / x; };
    double a2 = 1;
    double b2 = 2;
    double a_value2 = 1;
    double b_value2 = 0.5;
    model::Func precise_solution2 = [](double x) { return 1 / (x * x); };
    model::BoundaryCondition<model::BoundaryConditionKind::ThirdKind> cond4 = {a2, 1, 0, a_value2,
                                                                               b2, 3, 1, b_value2};
    model::LinearODESolver solver4({p2, q2}, f2, cond4);
    std::cout << "N = 1:" << '\n';
    auto collocation_approximate_solution6 = solver4.Solve(1, model::ODESolveMethod::Collocation);
    auto galerkin_approximate_solution6 = solver4.Solve(1, model::ODESolveMethod::Galerkin);
    semester6_task7::CompareSolutions(collocation_approximate_solution6,
                                      galerkin_approximate_solution6, precise_solution2, a2, b2);
    std::cout << "N = 5:" << '\n';
    auto collocation_approximate_solution7 = solver4.Solve(5, model::ODESolveMethod::Collocation);
    auto galerkin_approximate_solution7 = solver4.Solve(5, model::ODESolveMethod::Galerkin);
    semester6_task7::CompareSolutions(collocation_approximate_solution7,
                                      galerkin_approximate_solution7, precise_solution2, a2, b2);
    std::cout << "N = 10:" << '\n';
    auto collocation_approximate_solution8 = solver4.Solve(10, model::ODESolveMethod::Collocation);
    auto galerkin_approximate_solution8 = solver4.Solve(10, model::ODESolveMethod::Galerkin);
    semester6_task7::CompareSolutions(collocation_approximate_solution8,
                                      galerkin_approximate_solution8, precise_solution2, a2, b2);
    std::cout << "N = 20:" << '\n';
    auto collocation_approximate_solution9 = solver4.Solve(20, model::ODESolveMethod::Collocation);
    auto galerkin_approximate_solution9 = solver4.Solve(20, model::ODESolveMethod::Galerkin);
    semester6_task7::CompareSolutions(collocation_approximate_solution9,
                                      galerkin_approximate_solution9, precise_solution2, a2, b2);
}

}  // namespace tasks
