#pragma once

#include <cstddef>

#include "function.h"

namespace model {
class DefiniteIntegralCalculator {
private:
    Func f;

    double CalculateW(double a, double b, std::size_t n) const;
    double CalculateQ(double a, double b, std::size_t n) const;

public:
    DefiniteIntegralCalculator(Func fun) : f(fun) {}

    double Integral(double a, double b) const;
};

}  // namespace model
