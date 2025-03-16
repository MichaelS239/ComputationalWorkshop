#include "task4_3.h"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "calculate_integral.h"
#include "model/function.h"
#include "model/polynomial.h"
#include "util/input_util.h"

namespace semester5_task4_3 {
void PrintIntegrals(std::vector<std::pair<std::string, double>> const& integrals, double integral) {
    for (std::size_t i = 0; i != integrals.size(); ++i) {
        std::cout << "Approximate value of the intergal of the function obtained by the composite "
                  << integrals[i].first << ": " << std::setprecision(19) << integrals[i].second
                  << '\n';
        std::cout << "Absolute error: " << std::setprecision(19)
                  << std::abs(integrals[i].second - integral) << '\n';
        std::cout << "Relative error: " << std::setprecision(19)
                  << std::abs(integrals[i].second - integral) / std::abs(integral) << '\n';
    }
}

}  // namespace semester5_task4_3

namespace tasks {
void Semester5Task4_3() {
    std::cout << "Numerical integration using composite quadrature formulas" << '\n';
    std::cout << "Function to study: e^x * sin(x)" << '\n';
    std::cout << "Do you wish to find the integral of f(x) or of a random polynomial? [1|2]"
              << '\n';
    char c = util::InputChoice('1', '2');
    model::Polynomial polynomial;
    if (c == '2') {
        std::cout << "Enter the degree of the polynomial:" << '\n';
        int n;
        std::cin >> n;
        while (n < 0) {
            std::cout << "Degree must be non-negative." << '\n';
            std::cout << "Enter the degree of the polynomial:" << '\n';
            std::cin >> n;
        }
        n++;
        polynomial = model::Polynomial::CreateFromRandom(n);
        std::cout << "Polynomial: " << polynomial.ToString() << '\n';
    }
    auto [a, b] = util::InputBoundaries("integration segment");
    int n = util::InputSegments(a, b);

    model::Func f;
    double integral;
    if (c == '1') {
        f = semester5_task4_3::f;
        integral = semester5_task4_3::Integral(a, b);
    } else {
        f = polynomial;
        integral = polynomial.Integral(a, b);
    }
    std::cout << "Precise value of the intergal of the function: " << std::setprecision(15)
              << integral << '\n';
    std::cout << '\n';

    std::cout << "Approximate values for " << n << " partition segments:" << '\n';
    auto integrals = semester5_task4_3::CalculateIntegrals(a, b, n, f);
    semester5_task4_3::PrintIntegrals(integrals, integral);
    std::cout << '\n';

    std::cout << "Enter the number l in order to increase the number of partition segments l "
                 "times (l is a natural number):"
              << '\n';
    int l;
    std::cin >> l;
    while (l <= 1) {
        std::cout << "l must be greater than 1." << '\n';
        std::cout << "Enter the number l in order to increase the number of partition segments l "
                     "times (l is a natural number):"
                  << '\n';
        std::cin >> l;
    }
    n *= l;
    std::cout << "Approximate values for " << n << " partition segments:" << '\n';
    auto new_integrals = semester5_task4_3::CalculateIntegrals(a, b, n, f);
    semester5_task4_3::PrintIntegrals(new_integrals, integral);
    std::cout << '\n';
    std::cout << "Clarified values by the Runge rule:" << '\n';
    auto runge_integrals = semester5_task4_3::CalculateRunge(integrals, new_integrals, l);
    semester5_task4_3::PrintIntegrals(runge_integrals, integral);
}

}  // namespace tasks
