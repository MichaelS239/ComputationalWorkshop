#include "task6.h"

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
#include "model/root_calculator.h"
#include "util/create_system_matrix.h"
#include "util/matrix_util.h"
#include "util/table.h"

namespace tasks {
void Semester5Task6() {
    std::cout << "Numerical weighted integral calculation using quadrature formulas with maximum "
                 "order of accuracy"
              << '\n';
    std::cout << "Function to study: f(x) = sin(x)" << '\n';
    std::cout << "Segment: [0, 1]" << '\n';
    std::cout << "Precise value of the integral of the function: " << std::setprecision(15)
              << semester5_task6::integral << '\n';
    std::cout << '\n';
    std::cout << "Enter the number of points:" << '\n';
    int n;
    std::cin >> n;
    while (n <= 0) {
        std::cout << "The number of points must be positive." << '\n';
        std::cout << "Enter the number of points:" << '\n';
        std::cin >> n;
    }

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
    std::cout << "Weight function: ρ(x) = -xln(x)" << '\n';
    std::vector<double> moments = semester5_task6::CalculateMoments(n);
    std::vector<double> numbers = std::vector<double>(moments.size(), 0);
    std::iota(numbers.begin(), numbers.end(), 1);
    std::vector<std::pair<double, double>> moments_table = util::CreateTable(numbers, moments);
    std::cout << "Moments of the weight function:" << '\n';
    util::PrintTable(moments_table, {"k", "M_k"});
    std::vector<std::vector<double>> matrix = util::CreateSystemMatrix(points);
    std::vector<double> coefs = util::SolveSystem(matrix, moments);
    std::vector<std::pair<double, double>> table = util::CreateTable(points, coefs);
    std::cout << "The table of coefficients for every point:" << '\n';
    util::PrintTable(table, {"x_k", "A_k"});

    model::Polynomial polynomial = model::Polynomial::CreateFromRandom(n);
    double polynomial_integral = polynomial.WeightedIntegral(moments);
    std::cout << "Polynomial: " << polynomial.ToString() << '\n';
    model::IntegralCalculator polynomial_integral_calculator(table, polynomial);
    polynomial_integral_calculator.PrintDeviation(polynomial_integral);

    std::cout << "Function: f(x) = sin(x)" << '\n';
    model::IntegralCalculator function_integral_calculator(std::move(table), semester5_task6::f);
    function_integral_calculator.PrintDeviation(semester5_task6::integral);

    std::cout << "Quadrature formulas with maximum order of accuracy:" << '\n';
    std::cout << '\n';

    for (int weight = 0; weight != 3; ++weight) {
        double a;
        double b = 1;
        double integral;
        std::vector<double> old_moments;
        std::vector<double> new_moments;
        switch (weight) {
            case 0:
                std::cout << "Weight function: ρ(x) = -xln(x), segment [0, 1]" << '\n';
                a = 0;
                old_moments = semester5_task6::CalculateMoments(n);
                new_moments = semester5_task6::CalculateMoments(2 * n);
                integral = semester5_task6::integral;
                break;
            case 1:
                std::cout << "Weight function: ρ(x) = 1, segment [-1, 1]" << '\n';
                a = -1;
                old_moments = semester5_task6::CalculateGaussMoments(n);
                new_moments = semester5_task6::CalculateGaussMoments(2 * n);
                integral = 0;
                break;
            case 2:
                std::cout << "Weight function: ρ(x) = 1/√(1-x^2), segment [-1, 1]" << '\n';
                a = -1;
                old_moments = semester5_task6::CalculateMellerMoments(n);
                new_moments = semester5_task6::CalculateMellerMoments(2 * n);
                integral = 0;
                break;
            default:
                break;
        }

        numbers = std::vector<double>(new_moments.size(), 0);
        std::iota(numbers.begin(), numbers.end(), 1);
        std::vector<std::pair<double, double>> new_moments_table =
                util::CreateTable(numbers, new_moments);
        std::cout << "Moments of the weight function:" << '\n';
        util::PrintTable(new_moments_table, {"k", "M_k"});
        auto const [points_matrix, points_vector] =
                semester5_task6::CreatePointsSystem(n, new_moments);
        std::vector<double> solution = util::SolveSystem(points_matrix, points_vector);
        solution.push_back(1);
        model::Polynomial ortho_polynomial = std::move(solution);
        std::cout << "Orthogonal polynomial:" << '\n';
        std::cout << ortho_polynomial.ToString() << '\n';
        auto const& func = ortho_polynomial;
        auto const& func_prime = [ortho_polynomial](double x) { return ortho_polynomial.Prime(x); };
        model::RootCalculator root_calculator(func, func_prime);
        std::vector<double> new_points = root_calculator.FindRoots(a, b, 1000, 1e-8);
        numbers = std::vector<double>(new_points.size(), 0);
        std::iota(numbers.begin(), numbers.end(), 1);
        std::vector<std::pair<double, double>> points_table =
                util::CreateTable(numbers, new_points);
        std::cout << "Roots of the orthogonal polynomial:" << '\n';
        util::PrintTable(points_table, {"k", "x_k"});
        std::vector<std::vector<double>> coef_matrix = util::CreateSystemMatrix(new_points);
        std::vector<double> new_coefs = util::SolveSystem(coef_matrix, old_moments);
        std::vector<std::pair<double, double>> new_table = util::CreateTable(new_points, new_coefs);
        std::cout << "The table of coefficients for every point:" << '\n';
        util::PrintTable(new_table, {"x_k", "A_k"});

        model::Polynomial new_polynomial = model::Polynomial::CreateFromRandom(2 * n);
        double new_polynomial_integral = polynomial.WeightedIntegral(new_moments);
        std::cout << "Polynomial: " << new_polynomial.ToString() << '\n';
        model::IntegralCalculator new_polynomial_integral_calculator(new_table, new_polynomial);
        new_polynomial_integral_calculator.PrintDeviation(new_polynomial_integral);

        model::Polynomial polynomial2 = std::vector<double>(2 * n);
        polynomial2[0] = 1.125;
        polynomial2[1] = -2.55;
        polynomial2[polynomial2.Degree() - 1] = 0.175;
        double polynomial2_integral = polynomial2.WeightedIntegral(new_moments);
        std::cout << "Polynomial: " << polynomial2.ToString() << '\n';
        model::IntegralCalculator polynomial2_integral_calculator(new_table, polynomial2);
        polynomial2_integral_calculator.PrintDeviation(polynomial2_integral);

        std::cout << "Function: f(x) = sin(x)" << '\n';
        model::IntegralCalculator new_function_integral_calculator(std::move(new_table),
                                                                   semester5_task6::f);
        new_function_integral_calculator.PrintDeviation(integral);
    }
}

}  // namespace tasks
