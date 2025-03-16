#include "task1.h"

#include <iostream>
#include <utility>
#include <vector>

#include "calculate_roots.h"
#include "model/root_calculator.h"
#include "util/input_util.h"

namespace tasks {

void Semester5Task1() {
    std::cout << "Root-finding algorithms" << '\n';
    std::cout << "Function under study: f(x) = 8*cos(x)-x-6" << '\n';

    model::RootCalculator root_calculator(semester5_task1::f, semester5_task1::f_prime);

    auto [a, b] = util::InputBoundaries("root search");

    std::vector<std::pair<double, double>> root_segments;
    char c = 'y';
    while (c == 'y') {
        int n = util::InputSegments(a, b, 2);
        root_segments = root_calculator.FindRootSegments(a, b, n);
        std::cout << "There were found " << root_segments.size()
                  << " partition segments that contain roots of the equation f(x) = 0 on the "
                     "segment ["
                  << a << ", " << b << "]:" << '\n';
        for (auto const& [left, right] : root_segments) {
            std::cout << "[" << left << ", " << right << "]" << '\n';
        }
        std::cout << "Do you want to change the number of partition segments (if it increases, it "
                     "may be possible to find new roots)? [y|n]"
                  << '\n';
        c = util::InputChoice('y', 'n');
    }

    double eps;
    std::cout << "Enter the error bound of the root search (ε):" << '\n';
    std::cin >> eps;
    while (eps <= 0) {
        std::cout << "Error bound must be a positive value." << '\n';
        std::cout << "Enter the error bound of the root search (ε):" << '\n';
        std::cin >> eps;
    }

    semester5_task1::FindRoots(std::move(root_calculator), std::move(root_segments), eps);
}

}  // namespace tasks
