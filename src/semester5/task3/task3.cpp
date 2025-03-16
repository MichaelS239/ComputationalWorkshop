#include "task3.h"

#include <iostream>
#include <vector>

#include "calculate_derivatives.h"
#include "model/function.h"
#include "util/input_util.h"
#include "util/table.h"

namespace tasks {
void Semester5Task3() {
    std::cout << "Numerical differentiation" << '\n';
    std::cout << '\n';
    char c = 'c';
    while (c != 'q') {
        std::cout << "Choose a function to study [1|2]:" << '\n';
        std::cout << "1) f(x) = x^2 / (1 + x^2)" << '\n';
        std::cout << "2) f(x) = e^(3x)" << '\n';
        char num = util::InputChoice('1', '2');

        model::Func func;
        model::Func func_prime;
        model::Func func_second_prime;
        if (num == '1') {
            func = semester5_task3::f1;
            func_prime = semester5_task3::f1_prime;
            func_second_prime = semester5_task3::f1_prime2;
        } else {
            func = semester5_task3::f2;
            func_prime = semester5_task3::f2_prime;
            func_second_prime = semester5_task3::f2_prime2;
        }

        int m;
        std::cout << "Enter the number of known values (greater than or equal to 5):" << '\n';
        std::cin >> m;
        while (m <= 4) {
            std::cout << "The number of values must be greater than or equal to 5." << '\n';
            std::cout << "Enter the number of known values:" << '\n';
            std::cin >> m;
        }

        double x0;
        std::cout << "Enter the initial value:" << '\n';
        std::cin >> x0;
        double h = util::InputStep();

        std::vector<std::pair<double, double>> table =
                semester5_task3::CalculateTable(x0, h, m, func);
        std::cout << "Function table:" << '\n';
        util::PrintTable(table);
        std::cout << "Proceed to clarification using the Runge rule (r) or print the table of "
                     "derivatives (t):"
                  << '\n';
        c = util::InputChoice('r', 't');
        if (c == 'r') {
            auto const& [runge, second_runge] = semester5_task3::CalculateRunge(
                    table, x0, h, m, func, func_prime, func_second_prime);
            std::cout << "The table of first derivatives:" << '\n';
            util::PrintTable(runge, {"x", "f(x)", "J(h)", "|f'(x)-J(h)|", "J(h/2)",
                                     "|f'(x)-J(h/2)|", "J1(x)", "|f'(x)-J1(x)|"});
            std::cout << "The table of second derivatives:" << '\n';
            util::PrintTable(second_runge, {"x", "f(x)", "J(h)", "|f''(x)-J(h)|", "J(h/2)",
                                            "|f''(x)-J(h/2)|", "J1(x)", "|f''(x)-J1(x)|"});
        } else {
            std::vector<std::vector<double>> derivatives =
                    semester5_task3::CalculateDerivatives(table, h, func_prime, func_second_prime);
            std::cout << "The table of derivatives:" << '\n';
            util::PrintTable(derivatives, {"x", "f(x)", "f1'(x)", "|f'(x)-f1'(x)|", "f2'(x)",
                                           "|f'(x)-f2'(x)|", "f1''(x)", "|f''(x)-f1''(x)|"});
        }
        std::cout << "Choose a new function (c) or quit (q):" << '\n';
        c = util::InputChoice('c', 'q');
    }
}

}  // namespace tasks
