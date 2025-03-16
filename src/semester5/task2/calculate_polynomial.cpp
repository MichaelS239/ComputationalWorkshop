#include "calculate_polynomial.h"

#include <vector>

namespace semester5_task2 {
double f(double x) {
    return x * x / (1 + x * x);
}

std::vector<std::pair<double, double>> CalculateTable(double a, double b, int m) {
    std::vector<std::pair<double, double>> table;
    double h = (b - a) / (m - 1);
    for (int i = 0; i != m; ++i) {
        table.emplace_back(a, f(a));
        a += h;
    }
    return table;
}

double CalculatePolynomial(std::vector<std::pair<double, double>> const& table, double x, int n) {
    double value = 0;
    for (int i = 0; i != n + 1; ++i) {
        double cur_value = 1;
        for (int j = 0; j != n + 1; ++j) {
            if (j != i) {
                cur_value = cur_value * (x - table[j].first) / (table[i].first - table[j].first);
            }
        }
        value += cur_value * table[i].second;
    }
    return value;
}

}  // namespace semester5_task2
