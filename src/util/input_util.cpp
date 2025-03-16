#include "input_util.h"

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

namespace util {
std::vector<double> InputPoints(std::size_t n) {
    std::vector<double> points = std::vector<double>(n);
    std::cout << "Enter the points:" << '\n';
    for (std::size_t i = 0; i != n; ++i) {
        while (true) {
            std::cin >> points[i];
            bool flag = true;
            for (std::size_t j = 0; j != i; ++j) {
                if (points[i] == points[j]) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                break;
            } else {
                std::cout << "Points must not repeat. Enter the point again:" << '\n';
            }
        }
    }

    return points;
}

char InputChoice(char first_choice, char second_choice) {
    char c;
    std::cin >> c;
    while (c != first_choice && c != second_choice) {
        std::cout << "Enter '" << first_choice << "' or '" << second_choice << "'." << '\n';
        std::cin >> c;
    }
    return c;
}

int InputSegments(double a, double b, int minimum_number) {
    int n;
    std::cout << "Enter the number of partition segments for the segment [" << a << ", " << b
              << "]:" << '\n';
    std::cin >> n;
    while (n < minimum_number) {
        std::cout << "The number of partition segments must be greater than or equal to "
                  << minimum_number << "." << '\n';
        std::cout << "Enter the number of partition segments for the segment [" << a << ", " << b
                  << "]:" << '\n';
        std::cin >> n;
    }
    return n;
}

double InputStep() {
    double h;
    std::cout << "Enter the step:" << '\n';
    std::cin >> h;
    while (h <= 0) {
        std::cout << "Step must be a positive number." << '\n';
        std::cout << "Enter the step:" << '\n';
        std::cin >> h;
    }
    return h;
}

std::pair<double, double> InputBoundaries(std::string const& description) {
    double a, b;
    std::cout << "Enter " << description << " boundaries (A and B):" << '\n';
    std::cin >> a >> b;
    while (a > b) {
        std::cout << "The left boundary must be less than or equal to the right one." << '\n';
        std::cout << "Enter " << description << " boundaries (A and B):" << '\n';
        std::cin >> a >> b;
    }
    return {a, b};
}

}  // namespace util
