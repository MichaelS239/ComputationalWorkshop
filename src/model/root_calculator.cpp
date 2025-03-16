#include "root_calculator.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

namespace model {

std::vector<std::pair<double, double>> RootCalculator::FindRootSegments(double a, double b,
                                                                        int n) const {
    double h = (b - a) / n;
    std::vector<std::pair<double, double>> root_segments;
    for (int step_num = 0; step_num != n; ++step_num) {
        if (f(a) * f(a + h) <= 0) {
            root_segments.emplace_back(a, a + h);
        }
        a += h;
    }
    return root_segments;
}

double RootCalculator::Bisection(double a, double b, double eps, bool print_info) const {
    double c;
    unsigned count = 0;
    while (b - a > 2 * eps) {
        c = (a + b) / 2;
        if (f(a) * f(c) <= 0) {
            b = c;
        } else {
            a = c;
        }
        ++count;
    }
    double x = (a + b) / 2;
    if (print_info) {
        std::cout << "The root of the function: " << std::setprecision(15) << x << '\n';
        std::cout << "The number of iterations: " << count << '\n';
        std::cout << "The length of the last segment: " << std::setprecision(15) << b - a << '\n';
        std::cout << "Accuracy: " << std::setprecision(15) << std::abs(f(x)) << '\n';
    }
    return x;
}

double RootCalculator::Newton(double a, double b, double eps, bool print_info) const {
    unsigned count = 0;
    double x0 = (a + b) / 2;
    if (print_info) {
        std::cout << "Initial approximation: " << std::setprecision(15) << x0 << '\n';
    }
    double x1 = x0 - f(x0) / f_prime(x0);
    while (std::abs(x1 - x0) > eps) {
        x0 = x1;
        x1 = x1 - f(x1) / f_prime(x1);
        ++count;
    }
    if (print_info) {
        std::cout << "The root of the function: " << std::setprecision(15) << x1 << '\n';
        std::cout << "The number of iterations: " << count << '\n';
        std::cout << "The length of the last segment: " << std::setprecision(15)
                  << std::abs(x1 - x0) << '\n';
        std::cout << "Accuracy: " << std::setprecision(15) << std::abs(f(x1)) << '\n';
    }
    return x1;
}

double RootCalculator::ModifiedNewton(double a, double b, double eps, bool print_info) const {
    unsigned count = 0;
    double x0 = (a + b) / 2;
    if (print_info) {
        std::cout << "Initial approximation: " << std::setprecision(15) << x0 << '\n';
    }
    double x1 = x0;
    double x2 = x0 - f(x0) / f_prime(x0);
    while (std::abs(x2 - x1) > eps) {
        x1 = x2;
        x2 = x2 - f(x2) / f_prime(x0);
        ++count;
    }
    if (print_info) {
        std::cout << "The root of the function: " << std::setprecision(15) << x2 << '\n';
        std::cout << "The number of iterations: " << count << '\n';
        std::cout << "The length of the last segment: " << std::setprecision(15)
                  << std::abs(x2 - x1) << '\n';
        std::cout << "Accuracy: " << std::setprecision(15) << std::abs(f(x2)) << '\n';
    }
    return x2;
}

double RootCalculator::Secant(double a, double b, double eps, bool print_info) const {
    unsigned count = 0;
    double x0 = a + (b - a) / 3;
    double x1 = a + 2 * (b - a) / 3;
    if (print_info) {
        std::cout << "Initial approximation: " << std::setprecision(15) << x0 << " " << x1 << '\n';
    }
    double x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
    double x3;
    while (std::abs(x2 - x1) > eps) {
        x3 = x2 - f(x2) * (x2 - x1) / (f(x2) - f(x1));
        x0 = x1;
        x1 = x2;
        x2 = x3;
        ++count;
    }
    if (print_info) {
        std::cout << "The root of the function: " << std::setprecision(15) << x2 << '\n';
        std::cout << "The number of iterations: " << count << '\n';
        std::cout << "The length of the last segment: " << std::setprecision(15)
                  << std::abs(x2 - x1) << '\n';
        std::cout << "Accuracy: " << std::setprecision(15) << std::abs(f(x2)) << '\n';
    }
    return x2;
}

std::vector<double> RootCalculator::FindRoots(double a, double b, double n, double eps) const {
    std::vector<std::pair<double, double>> root_segments = FindRootSegments(a, b, n);
    std::vector<double> roots;
    roots.reserve(root_segments.size());
    for (auto const& [left, right] : root_segments) {
        roots.push_back(Newton(left, right, eps));
    }
    return roots;
}

}  // namespace model
