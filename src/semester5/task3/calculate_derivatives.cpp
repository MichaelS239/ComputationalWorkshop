#include "calculate_derivatives.h"

#include <cmath>
#include <cstddef>
#include <vector>

namespace semester5_task3 {
double f1(double x) {
    return x * x / (1 + x * x);
}

double f1_prime(double x) {
    return (2 * x) / ((1 + x * x) * (1 + x * x));
}

double f1_prime2(double x) {
    return (2 * (1 + x * x) * (1 - 3 * x * x)) / (std::pow(1 + x * x, 4));
}

double f2(double x) {
    return std::exp(3 * x);
}

double f2_prime(double x) {
    return 3 * std::exp(3 * x);
}

double f2_prime2(double x) {
    return 9 * std::exp(3 * x);
}

std::vector<std::pair<double, double>> CalculateTable(double x0, double h, int m,
                                                      model::Func func) {
    std::vector<std::pair<double, double>> table;
    for (std::size_t i = 0; i != m; ++i) {
        table.emplace_back(x0, func(x0));
        x0 += h;
    }
    return table;
}

std::vector<std::vector<double>> CalculateFirstDerivatives(
        std::vector<std::pair<double, double>> const& table, double h, model::Func func_prime) {
    std::vector<std::vector<double>> derivatives;
    for (std::size_t i = 0; i != table.size(); ++i) {
        derivatives.emplace_back();
        derivatives[derivatives.size() - 1] = {0, 0};
        double derivative = 0;
        if (i == 0) {
            derivative = (-3 * table[i].second + 4 * table[i + 1].second - table[i + 2].second) /
                         (2 * h);
        } else if (i == table.size() - 1) {
            derivative =
                    (3 * table[i].second - 4 * table[i - 1].second + table[i - 2].second) / (2 * h);
        } else {
            derivative = (table[i + 1].second - table[i - 1].second) / (2 * h);
        }
        derivatives[derivatives.size() - 1][0] = derivative;
        derivatives[derivatives.size() - 1][1] = std::abs(func_prime(table[i].first) - derivative);
    }
    return derivatives;
}

std::vector<std::vector<double>> CalculateNewFirstDerivatives(
        std::vector<std::pair<double, double>> const& table, double h, model::Func func_prime) {
    std::vector<std::vector<double>> derivatives;
    for (std::size_t i = 0; i != table.size(); ++i) {
        derivatives.emplace_back();
        derivatives[derivatives.size() - 1] = {0, 0};
        double derivative = 0;
        if (i == 0) {
            derivative =
                    (-25 * table[i].second + 48 * table[i + 1].second - 36 * table[i + 2].second +
                     16 * table[i + 3].second - 3 * table[i + 4].second) /
                    (12 * h);
        } else if (i == 1) {
            derivative =
                    (-3 * table[i - 1].second - 10 * table[i].second + 18 * table[i + 1].second -
                     6 * table[i + 2].second + table[i + 3].second) /
                    (12 * h);
        } else if (i == table.size() - 2) {
            derivative =
                    (3 * table[i + 1].second + 10 * table[i].second - 18 * table[i - 1].second +
                     6 * table[i - 2].second - table[i - 3].second) /
                    (12 * h);
        } else if (i == table.size() - 1) {
            derivative =
                    (25 * table[i].second - 48 * table[i - 1].second + 36 * table[i - 2].second -
                     16 * table[i - 3].second + 3 * table[i - 4].second) /
                    (12 * h);
        } else {
            derivative = (table[i - 2].second - 8 * table[i - 1].second + 8 * table[i + 1].second -
                          table[i + 2].second) /
                         (12 * h);
        }
        derivatives[derivatives.size() - 1][0] = derivative;
        derivatives[derivatives.size() - 1][1] = std::abs(func_prime(table[i].first) - derivative);
    }
    return derivatives;
}

std::vector<std::vector<double>> CalculateSecondDerivatives(
        std::vector<std::pair<double, double>> const& table, double h,
        model::Func func_second_prime) {
    std::vector<std::vector<double>> derivatives;
    for (std::size_t i = 0; i != table.size(); ++i) {
        derivatives.emplace_back();
        derivatives[derivatives.size() - 1] = {0, 0};
        double derivative = 0;
        if (i == 0) {
            derivative = (2 * table[i].second - 5 * table[i + 1].second + 4 * table[i + 2].second -
                          table[i + 3].second) /
                         (h * h);
        } else if (i == table.size() - 1) {
            derivative = (2 * table[i].second - 5 * table[i - 1].second + 4 * table[i - 2].second -
                          table[i - 3].second) /
                         (h * h);
        } else {
            derivative =
                    (table[i + 1].second - 2 * table[i].second + table[i - 1].second) / (h * h);
        }
        derivatives[derivatives.size() - 1][0] = derivative;
        derivatives[derivatives.size() - 1][1] =
                std::abs(func_second_prime(table[i].first) - derivative);
    }
    return derivatives;
}

std::vector<std::vector<double>> CalculateDerivatives(
        std::vector<std::pair<double, double>> const& table, double h, model::Func func_prime,
        model::Func func_second_prime) {
    std::vector<std::vector<double>> derivatives = CalculateFirstDerivatives(table, h, func_prime);
    std::vector<std::vector<double>> new_derivatives =
            CalculateNewFirstDerivatives(table, h, func_prime);
    std::vector<std::vector<double>> second_derivatives =
            CalculateSecondDerivatives(table, h, func_second_prime);
    std::vector<std::vector<double>> derivatives_table(table.size());
    for (std::size_t i = 0; i != derivatives.size(); ++i) {
        derivatives_table[i].push_back(table[i].first);
        derivatives_table[i].push_back(table[i].second);
        for (std::size_t j = 0; j != derivatives[i].size(); ++j) {
            derivatives_table[i].push_back(derivatives[i][j]);
        }
        for (std::size_t j = 0; j != new_derivatives[i].size(); ++j) {
            derivatives_table[i].push_back(new_derivatives[i][j]);
        }
        for (std::size_t j = 0; j != second_derivatives[i].size(); ++j) {
            derivatives_table[i].push_back(second_derivatives[i][j]);
        }
    }
    return derivatives_table;
}

std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> CalculateRunge(
        std::vector<std::pair<double, double>> const& table, double x0, double h, int m,
        model::Func func, model::Func func_prime, model::Func func_second_prime) {
    std::vector<std::vector<double>> derivatives = CalculateFirstDerivatives(table, h, func_prime);
    std::vector<std::pair<double, double>> new_table = CalculateTable(x0, h / 2, 2 * m, func);
    std::vector<std::vector<double>> other_derivatives =
            CalculateFirstDerivatives(new_table, h / 2, func_prime);
    std::vector<std::vector<double>> derivatives_table(table.size());
    std::vector<double> runge(table.size());
    for (std::size_t i = 0; i != derivatives.size(); ++i) {
        runge[i] = (4 * other_derivatives[2 * i][0] - derivatives[i][0]) / 3;
    }
    for (std::size_t i = 0; i != derivatives.size(); ++i) {
        derivatives_table[i].push_back(table[i].first);
        derivatives_table[i].push_back(table[i].second);
        for (std::size_t j = 0; j != derivatives[i].size(); ++j) {
            derivatives_table[i].push_back(derivatives[i][j]);
        }
        for (std::size_t j = 0; j != other_derivatives[i].size(); ++j) {
            derivatives_table[i].push_back(other_derivatives[2 * i][j]);
        }
        derivatives_table[i].push_back(runge[i]);
        derivatives_table[i].push_back(std::abs(runge[i] - func_prime(table[i].first)));
    }

    std::vector<std::vector<double>> second_derivatives =
            CalculateSecondDerivatives(table, h, func_second_prime);
    std::vector<std::vector<double>> other_second_derivatives =
            CalculateSecondDerivatives(new_table, h / 2, func_second_prime);
    std::vector<std::vector<double>> second_derivatives_table(table.size());
    std::vector<double> second_runge(table.size());
    for (std::size_t i = 0; i != second_derivatives.size(); ++i) {
        second_runge[i] = (4 * other_second_derivatives[2 * i][0] - second_derivatives[i][0]) / 3;
    }
    for (std::size_t i = 0; i != second_derivatives.size(); ++i) {
        second_derivatives_table[i].push_back(table[i].first);
        second_derivatives_table[i].push_back(table[i].second);
        for (std::size_t j = 0; j != second_derivatives[i].size(); ++j) {
            second_derivatives_table[i].push_back(second_derivatives[i][j]);
        }
        for (std::size_t j = 0; j != other_second_derivatives[i].size(); ++j) {
            second_derivatives_table[i].push_back(other_second_derivatives[2 * i][j]);
        }
        second_derivatives_table[i].push_back(runge[i]);
        second_derivatives_table[i].push_back(
                std::abs(second_runge[i] - func_second_prime(table[i].first)));
    }
    return {derivatives_table, second_derivatives_table};
}

}  // namespace semester5_task3
