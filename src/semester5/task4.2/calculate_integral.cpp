#include "calculate_integral.h"

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace semester5_task4_2 {

double Antiderivative(double x) {
    return std::exp(x) * (std::sin(x) - std::cos(x)) / 2;
}

double Integral(double a, double b) {
    return Antiderivative(b) - Antiderivative(a);
}

double f(double x) {
    return std::exp(x) * std::sin(x);
}

double LeftRectangle(double a, double b, model::Func f) {
    return (b - a) * f(a);
}

double RightRectangle(double a, double b, model::Func f) {
    return (b - a) * f(b);
}

double MidpointRectangle(double a, double b, model::Func f) {
    return (b - a) * f((a + b) / 2);
}

double Trapezium(double a, double b, model::Func f) {
    return (b - a) / 2 * (f(a) + f(b));
}

double Simpson(double a, double b, model::Func f) {
    return (b - a) / 6 * (f(a) + 4 * f((a + b) / 2) + f(b));
}

double ThreeEigths(double a, double b, model::Func f) {
    return (b - a) / 8 * (f(a) + 3 * f(a + (b - a) / 3) + 3 * f(a + 2 * (b - a) / 3) + f(b));
}

std::vector<std::pair<std::string, double>> CalculateIntegrals(double a, double b, model::Func f) {
    std::vector<std::pair<std::string, double>> integrals;
    integrals.emplace_back("left rectangle rule", LeftRectangle(a, b, f));
    integrals.emplace_back("right rectangle rule", RightRectangle(a, b, f));
    integrals.emplace_back("midpoint rectangle rule", MidpointRectangle(a, b, f));
    integrals.emplace_back("trapezium rule", Trapezium(a, b, f));
    integrals.emplace_back("Simpson's 1/3 rule", Simpson(a, b, f));
    integrals.emplace_back("Simpson's 3/8 rule", ThreeEigths(a, b, f));
    return integrals;
}

}  // namespace semester5_task4_2
