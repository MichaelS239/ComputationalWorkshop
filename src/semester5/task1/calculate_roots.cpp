#include "calculate_roots.h"

#include <cmath>
#include <iostream>
#include <vector>

#include "model/root_calculator.h"

namespace semester5_task1 {
double f(double x) {
    return 8 * cos(x) - x - 6;
}

double f_prime(double x) {
    return -8 * sin(x) - 1;
}

void FindRoots(model::RootCalculator root_calculator,
               std::vector<std::pair<double, double>> root_segments, double eps) {
    for (auto const& [left, right] : root_segments) {
        std::cout << "Segment [" << left << ", " << right << "]:" << '\n';
        std::cout << '\n';
        std::cout << "Bisection method:" << '\n';
        root_calculator.Bisection(left, right, eps, true);
        std::cout << '\n';
        std::cout << "Newton's method:" << '\n';
        root_calculator.Newton(left, right, eps, true);
        std::cout << '\n';
        std::cout << "Modified Newton's method:" << '\n';
        root_calculator.ModifiedNewton(left, right, eps, true);
        std::cout << '\n';
        std::cout << "Secant method:" << '\n';
        root_calculator.Secant(left, right, eps, true);
        std::cout << '\n';
    }
}

}  // namespace semester5_task1
