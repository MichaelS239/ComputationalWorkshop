#include "definite_integral_calculator.h"

#include <cstddef>

namespace model {

double DefiniteIntegralCalculator::CalculateW(double a, double b, std::size_t n) const {
    double sum = 0;
    double h = (b - a) / n;
    a += h;
    for (std::size_t i = 1; i != n; ++i) {
        sum += f(a);
        a += h;
    }
    return sum;
}

double DefiniteIntegralCalculator::CalculateQ(double a, double b, std::size_t n) const {
    double sum = 0;
    double h = (b - a) / n;
    a += h / 2;
    for (std::size_t i = 0; i != n; ++i) {
        sum += f(a);
        a += h;
    }
    return sum;
}

double DefiniteIntegralCalculator::Integral(double a, double b) const {
    std::size_t n = 100;
    double h = (b - a) / n;
    double w = CalculateW(a, b, n);
    double q = CalculateQ(a, b, n);

    return h / 6 * (f(a) + f(b) + 2 * w + 4 * q);
}

}  // namespace model
