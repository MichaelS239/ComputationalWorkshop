#include "task5.h"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "calculate_equation.h"
#include "util/input_util.h"
#include "util/table.h"

namespace tasks {
void Semester5Task5() {
    std::cout << "Numerical methods for solving initial value problems for ODE" << '\n';
    std::cout << '\n';
    std::cout << "Initial value problem to study: y'(x) = -y(x) + e^(-x), y(0) = 1" << '\n';

    char c = 'y';
    while (c == 'y') {
        std::cout << "Enter the number of points N (N >= 3):" << '\n';
        int n;
        std::cin >> n;
        while (n < 3) {
            std::cout << "N must be greater than or equal to 3." << '\n';
            std::cout << "Enter N:" << '\n';
            std::cin >> n;
        }
        double h = util::InputStep();

        std::vector<std::vector<double>> taylor_table = semester5_task5::CompareWithSolution(
                -2, n + 1, h, semester5_task5::CalculateTaylor);
        std::cout << "The table of results for Taylor series method:" << '\n';
        util::PrintTable(taylor_table, {"x", "y(x)", "T(x)", "|y(x)-T(x)|"});
        std::cout << '\n';

        std::vector<std::vector<double>> adams_table =
                semester5_task5::CompareWithSolution(3, n + 1, h, semester5_task5::CalculateAdams);
        std::cout << "The table of results for Adams' method:" << '\n';
        util::PrintTable(adams_table, {"x", "y(x)", "A(x)", "|y(x)-A(x)|"});
        std::cout << '\n';

        std::vector<std::vector<double>> runge_kutta_table = semester5_task5::CompareWithSolution(
                1, n + 1, h, semester5_task5::CalculateRungeKutta);
        std::cout << "The table of results for Runge-Kutta method:" << '\n';
        util::PrintTable(runge_kutta_table, {"x", "y(x)", "R(x)", "|y(x)-R(x)|"});
        std::cout << '\n';

        std::vector<std::vector<double>> euler_table =
                semester5_task5::CompareWithSolution(1, n + 1, h, semester5_task5::CalculateEuler);
        std::cout << "The table of results for Euler's method:" << '\n';
        util::PrintTable(euler_table, {"x", "y(x)", "E(x)", "|y(x)-E(x)|"});
        std::cout << '\n';

        std::vector<std::vector<double>> euler1_table =
                semester5_task5::CompareWithSolution(1, n + 1, h, semester5_task5::CalculateEuler1);
        std::cout << "The table of results for Euler's method 1:" << '\n';
        util::PrintTable(euler1_table, {"x", "y(x)", "E1(x)", "|y(x)-E1(x)|"});
        std::cout << '\n';

        std::vector<std::vector<double>> euler2_table =
                semester5_task5::CompareWithSolution(1, n + 1, h, semester5_task5::CalculateEuler2);
        std::cout << "The table of results for Euler's method 2:" << '\n';
        util::PrintTable(euler2_table, {"x", "y(x)", "E2(x)", "|y(x)-E2(x)|"});
        std::cout << '\n';

        std::cout << "The table of absolute errors for the last point:" << '\n';
        std::vector<std::vector<double>> final_table =
                std::vector<std::vector<double>>(1, std::vector<double>(8));
        final_table[0][0] = n * h;
        final_table[0][1] = semester5_task5::Solution(n * h);
        final_table[0][2] = taylor_table[taylor_table.size() - 1][3];
        final_table[0][3] = adams_table[adams_table.size() - 1][3];
        final_table[0][4] = runge_kutta_table[runge_kutta_table.size() - 1][3];
        final_table[0][5] = euler_table[euler_table.size() - 1][3];
        final_table[0][6] = euler1_table[euler1_table.size() - 1][3];
        final_table[0][7] = euler2_table[euler2_table.size() - 1][3];
        util::PrintTable(final_table, {"x", "y(x)", "|y(x)-T(x)|", "|y(x)-A(x)|", "|y(x)-R(x)|",
                                       "|y(x)-E(x)|", "|y(x)-E1(x)|", "|y(x)-E2(x)|"});
        std::cout << '\n';

        std::cout << "Do you want to enter a new number of points and a new step? [y|n]" << '\n';
        c = util::InputChoice('y', 'n');
    }
}

}  // namespace tasks
