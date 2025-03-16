#include "task1.h"

#include <iostream>
#include <utility>
#include <vector>

#include "calculate_roots.h"
#include "model/root_calculator.h"

namespace semester5_task1 {
int InputSegments(double a, double b) {
    int n;
    std::cout << "Enter the number of partition segments for the segment [" << a << ", " << b
              << "]:" << '\n';
    std::cin >> n;
    while (n < 2) {
        std::cout << "The number of partition segments must be greater than or equal to two."
                  << '\n';
        std::cout << "Enter the number of partition segments for the segment [" << a << ", " << b
                  << "]:" << '\n';
        std::cin >> n;
    }
    return n;
}

}  // namespace semester5_task1

namespace tasks {

void Semester5Task1() {
    std::cout << "Root-finding algorithms" << '\n';
    std::cout << "Function under study: f(x) = 8*cos(x)-x-6" << '\n';

    model::RootCalculator root_calculator(semester5_task1::f, semester5_task1::f_prime);

    double a, b;
    std::cout << "Enter root search boundaries (A and B):" << '\n';
    std::cin >> a >> b;
    while (a > b) {
        std::cout << "The left boundary must be less than or equal to the right one." << '\n';
        std::cout << "Enter root search boundaries (A and B):" << '\n';
        std::cin >> a >> b;
    }

    std::vector<std::pair<double, double>> root_segments;
    char c = 'y';
    while (c == 'y') {
        int n = semester5_task1::InputSegments(a, b);
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
        std::cin >> c;
        while (c != 'y' && c != 'n') {
            std::cout << "Enter 'y' or 'n'." << '\n';
            std::cin >> c;
        }
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
