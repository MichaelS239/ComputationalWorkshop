#include "task4_2.h"

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <vector>

#include "calculate_integral.h"
#include "model/function.h"
#include "model/polynomial.h"
#include "util/input_util.h"

namespace tasks {
void Semester5Task4_2() {
    std::cout << "Numerical integration using simple quadrature formulas" << '\n';
    std::cout << '\n';
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

    model::Func f;
    double integral;
    if (c == '1') {
        f = semester5_task4_2::f;
        integral = semester5_task4_2::Integral(a, b);
    } else {
        f = polynomial;
        integral = polynomial.Integral(a, b);
    }
    std::cout << "Precise value of the intergal of the function: " << std::setprecision(15)
              << integral << '\n';
    std::cout << '\n';

    auto integrals = semester5_task4_2::CalculateIntegrals(a, b, f);
    for (std::size_t i = 0; i != integrals.size(); ++i) {
        std::cout << "Approximate value of the intergal of the function obtained by the "
                  << integrals[i].first << ": " << std::setprecision(15) << integrals[i].second
                  << '\n';
        std::cout << "Absolute error: " << std::abs(integrals[i].second - integral) << '\n';
    }
}

}  // namespace tasks
