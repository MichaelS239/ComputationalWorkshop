#include "task2.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "calculate_polynomial.h"
#include "util/input_util.h"
#include "util/table.h"

namespace tasks {
void Semester5Task2() {
    std::cout << "Polynomial interpolation" << '\n';
    std::cout << "Function under study: x^2 / (1 + x^2)" << '\n';

    int m;
    std::cout << "Enter the number of known values:" << '\n';
    std::cin >> m;
    while (m <= 0) {
        std::cout << "The number of known values must be positive." << '\n';
        std::cout << "Enter the number of known values:" << '\n';
        std::cin >> m;
    }

    auto [a, b] = util::InputBoundaries("interpolation segment");

    std::vector<std::pair<double, double>> table = semester5_task2::CalculateTable(a, b, m);
    std::cout << "The table of function values:" << '\n';
    util::PrintTable(table);

    char c = 'y';
    while (c == 'y') {
        double x;
        std::cout << "Enter the point of interpolation (x), at which you wish to find the value of "
                     "the function:"
                  << '\n';
        std::cin >> x;

        int n;
        std::cout << "Enter the degree of the interpolating polynomial n (n <= " << m - 1
                  << "):" << '\n';
        std::cin >> n;
        while (n <= 0 || n >= m) {
            std::cout << "Invalid value entered." << '\n';
            std::cout << "Enter the degree of the interpolating polynomial n (n <= " << m - 1
                      << "):" << '\n';
            std::cin >> n;
        }

        std::sort(table.begin(), table.end(), [x](auto const& left_pair, auto const& right_pair) {
            return std::abs(left_pair.first - x) < std::abs(right_pair.first - x);
        });
        std::cout << "Sorted function table:" << '\n';
        util::PrintTable(table);

        double value = semester5_task2::CalculatePolynomial(table, x, n);
        std::cout << "The value of the Lagrange interpolating polynomial at the point "
                  << std::setprecision(15) << x << ": " << value << '\n';
        double dif = std::abs(semester5_task2::f(x) - value);
        std::cout << "Error of interpolation: " << std::setprecision(15) << dif << '\n';
        std::cout << '\n';
        std::cout << "Do you want to enter new values of x and n? [y|n]" << '\n';
        c = util::InputChoice('y', 'n');
    }
}

}  // namespace tasks
