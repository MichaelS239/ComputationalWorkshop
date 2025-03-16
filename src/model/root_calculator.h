#pragma once

#include <vector>

#include "function.h"

namespace model {
class RootCalculator {
private:
    Func f;
    Func f_prime;

public:
    RootCalculator(Func func, Func func_prime) : f(func), f_prime(func_prime) {}

    std::vector<std::pair<double, double>> FindRootSegments(double a, double b, int n) const;

    double Bisection(double a, double b, double eps, bool print_info = false) const;
    double Newton(double a, double b, double eps, bool print_info = false) const;
    double ModifiedNewton(double a, double b, double eps, bool print_info = false) const;
    double Secant(double a, double b, double eps, bool print_info = false) const;

    std::vector<double> FindRoots(double a, double b, double n, double eps) const;
};

}  // namespace model
