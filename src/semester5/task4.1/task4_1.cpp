#include "task4_1.h"

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "calculate_integral.h"
#include "model/integral_calculator.h"
#include "model/polynomial.h"
#include "util/input_util.h"
#include "util/matrix_util.h"
#include "util/table.h"

namespace tasks {
void Semester5Task4_1() {
    std::cout << "Numerical weighted integral calculation" << '\n';
    std::cout << '\n';
    std::cout << "Function to study: f(x) = sin(x)" << '\n';
    std::cout << "Segment: [0, 1]" << '\n';
    std::cout << "Precise value of the integral: " << std::setprecision(15)
              << semester5_task4_1::integral << '\n';
    std::cout << '\n';
    std::cout << "Enter the number of points:" << '\n';
    int n;
    std::cin >> n;
    while (n <= 0) {
        std::cout << "The number of points must be positive." << '\n';
        std::cout << "Enter the number of points:" << '\n';
        std::cin >> n;
    }

    std::vector<double> points = util::InputPoints(n);
    std::cout << "Weight function: ρ(x) = -xln(x)" << '\n';
    auto [matrix, vector] = semester5_task4_1::CreateSystem(points);
    std::vector<double> numbers = std::vector<double>(vector.size(), 0);
    std::iota(numbers.begin(), numbers.end(), 1);
    std::vector<std::pair<double, double>> moments_table = util::CreateTable(numbers, vector);
    std::cout << "Moments of the weight function:" << '\n';
    util::PrintTable(moments_table, {"k", "M_k"});
    std::vector<double> coefs = util::SolveSystem(matrix, vector);
    std::vector<std::pair<double, double>> table = util::CreateTable(points, coefs);
    std::cout << "The table of coefficients for every point:" << '\n';
    util::PrintTable(table, {"x_k", "A_k"});

    model::Polynomial polynomial = model::Polynomial::CreateFromRandom(n);
    std::vector<double> moments = std::vector<double>(n);
    for (std::size_t i = 0; i != n; ++i) {
        moments[i] = 1.0 / ((i + 2) * (i + 2));
    }
    double polynomial_integral = polynomial.WeightedIntegral(moments);
    std::cout << "Polynomial: " << polynomial.ToString() << '\n';
    model::IntegralCalculator polynomial_integral_calculator(table, polynomial);
    polynomial_integral_calculator.PrintDeviation(polynomial_integral);

    std::cout << "Function: f(x) = sin(x)" << '\n';
    model::IntegralCalculator function_integral_calculator(std::move(table), semester5_task4_1::f);
    function_integral_calculator.PrintDeviation(semester5_task4_1::integral);
}

}  // namespace tasks
