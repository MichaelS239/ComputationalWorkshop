#pragma once

#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

#include "function.h"
#include "util/table.h"

namespace model {
class IntegralCalculator {
private:
    std::vector<std::pair<double, double>> values;
    Func f;

public:
    IntegralCalculator(std::vector<std::pair<double, double>> const& table, Func func)
        : values(table), f(func) {}

    IntegralCalculator(std::vector<std::pair<double, double>>&& table, Func func)
        : values(std::move(table)), f(func) {}

    IntegralCalculator(std::vector<double> const& points, std::vector<double> const& coefs,
                       Func func)
        : values(util::CreateTable(points, coefs)), f(func) {}

    IntegralCalculator(std::vector<double>&& points, std::vector<double>&& coefs, Func func)
        : values(util::CreateTable(std::move(points), std::move(coefs))), f(func) {}

    double CalculateIntegral() const {
        double sum = 0;
        for (std::size_t i = 0; i != values.size(); ++i) {
            sum += values[i].second * f(values[i].first);
        }
        return sum;
    }

    void PrintDeviation(double precise_integral) {
        std::cout << "Precise value of the intergal of the function: " << std::setprecision(15)
                  << precise_integral << '\n';
        double approximate_integral = CalculateIntegral();
        std::cout << "Approximate value of the intergal of the function: " << std::setprecision(15)
                  << approximate_integral << '\n';
        std::cout << "Absolute error: " << std::setprecision(15)
                  << std::abs(approximate_integral - precise_integral) << '\n';
        std::cout << '\n';
    }
};
}  // namespace model
