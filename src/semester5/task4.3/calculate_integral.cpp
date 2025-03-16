#include "calculate_integral.h"

#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace semester5_task4_3 {

double Antiderivative(double x) {
    return std::exp(x) * (std::sin(x) - std::cos(x)) / 2;
}

double Integral(double a, double b) {
    return Antiderivative(b) - Antiderivative(a);
}

double f(double x) {
    return std::exp(x) * std::sin(x);
}

double CalculateW(double a, double b, int n, model::Func f) {
    double sum = 0;
    double h = (b - a) / n;
    a += h;
    for (int i = 1; i != n; ++i) {
        sum += f(a);
        a += h;
    }
    return sum;
}

double LeftRectangle(double a, double w, double h, model::Func f) {
    return h * (f(a) + w);
}

double RightRectangle(double b, double w, double h, model::Func f) {
    return h * (w + f(b));
}

double CalculateQ(double a, double b, int n, model::Func f) {
    double sum = 0;
    double h = (b - a) / n;
    a += h / 2;
    for (int i = 0; i != n; ++i) {
        sum += f(a);
        a += h;
    }
    return sum;
}

double MidpointRectangle(double h, double q) {
    return h * q;
}

double Trapezium(double a, double b, double w, double h, model::Func f) {
    return h / 2 * (f(a) + f(b) + 2 * w);
}

double Simpson(double a, double b, double w, double q, double h, model::Func f) {
    return h / 6 * (f(a) + f(b) + 2 * w + 4 * q);
}

std::vector<std::pair<std::string, double>> CalculateIntegrals(double a, double b, int n,
                                                               model::Func f) {
    std::vector<std::pair<std::string, double>> integrals;
    double h = (b - a) / n;
    double w = CalculateW(a, b, n, f);
    double q = CalculateQ(a, b, n, f);
    integrals.emplace_back("left rectangle rule", LeftRectangle(a, w, h, f));
    integrals.emplace_back("right rectangle rule", RightRectangle(b, w, h, f));
    integrals.emplace_back("midpoint rectangle rule", MidpointRectangle(h, q));
    integrals.emplace_back("trapezium rule", Trapezium(a, b, w, h, f));
    integrals.emplace_back("Simpson's 1/3 rule", Simpson(a, b, w, q, h, f));
    return integrals;
}

std::vector<std::pair<std::string, double>> CalculateRunge(
        std::vector<std::pair<std::string, double>> const& integrals,
        std::vector<std::pair<std::string, double>> const& new_integrals, int l) {
    std::vector<std::pair<std::string, double>> runge_integrals;
    for (std::size_t i = 0; i != integrals.size(); ++i) {
        runge_integrals.emplace_back(
                integrals[i].first,
                (std::pow(l, precisions[i]) * new_integrals[i].second - integrals[i].second) /
                        (std::pow(l, precisions[i]) - 1));
    }
    return runge_integrals;
}

}  // namespace semester5_task4_3
