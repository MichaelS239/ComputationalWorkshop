#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "boundary_condition.h"
#include "definite_integral_calculator.h"
#include "function.h"
#include "matrix.h"
#include "solve_methods.h"

namespace model {

using SolveInfo = std::vector<std::pair<double, double>>;

template <BoundaryConditionKind T>
class LinearODESolver {
private:
    std::pair<Func, Func> lhs_;
    Func rhs_;
    BoundaryCondition<T> boundary_condition_;

    std::vector<double> CreatePoints(double a, double b, double h, std::size_t n) const {
        std::vector<double> points(n + 1);
        for (std::size_t i = 0; i != n + 1; ++i) {
            points[i] = a + i * h;
        }
        return points;
    }

    std::vector<double> CreateInnerPoints(double a, double b, std::size_t n) const {
        std::vector<double> points(n);
        double h = (b - a) / (n + 1);
        for (std::size_t i = 0; i != n; ++i) {
            points[i] = a + (i + 1) * h;
        }
        return points;
    }

    std::vector<double> CalculateSolution(std::size_t points_num) const {
        double a = boundary_condition_.left_boundary;
        double b = boundary_condition_.right_boundary;
        double h = (b - a) / points_num;
        std::vector<double> points = CreatePoints(a, b, h, points_num);
        Matrix matrix(points_num + 1);
        std::vector<double> coefs(points_num + 1);
        if constexpr (T == BoundaryConditionKind::FirstKind) {
            matrix[0][0] = 1;
            coefs[0] = boundary_condition_.left_value;
            matrix[points_num][points_num] = 1;
            coefs[points_num] = boundary_condition_.right_value;
        } else if constexpr (T == BoundaryConditionKind::SecondKind) {
            matrix[0][0] = -3 / (2 * h);
            matrix[0][1] = 2 / h;
            matrix[0][2] = -1 / (2 * h);
            coefs[0] = boundary_condition_.left_value;
            matrix[points_num][points_num] = 3 / (2 * h);
            matrix[points_num][points_num - 1] = -2 / h;
            matrix[points_num][points_num - 2] = 1 / (2 * h);
            coefs[points_num] = boundary_condition_.right_value;
        } else if constexpr (T == BoundaryConditionKind::ThirdKind) {
            matrix[0][0] = boundary_condition_.left_value_coef -
                           3 * boundary_condition_.left_derivative_coef / (2 * h);
            matrix[0][1] = 2 * boundary_condition_.left_derivative_coef / h;
            matrix[0][2] = -boundary_condition_.left_derivative_coef / (2 * h);
            coefs[0] = boundary_condition_.left_value;
            matrix[points_num][points_num] =
                    boundary_condition_.right_value_coef +
                    3 * boundary_condition_.right_derivative_coef / (2 * h);
            matrix[points_num][points_num - 1] = -2 * boundary_condition_.right_derivative_coef / h;
            matrix[points_num][points_num - 2] =
                    boundary_condition_.right_derivative_coef / (2 * h);
            coefs[points_num] = boundary_condition_.right_value;
        }
        for (std::size_t i = 1; i != points_num; ++i) {
            matrix[i][i - 1] = 1 / (h * h) - lhs_.first(points[i]) / (2 * h);
            matrix[i][i] = -2 / (h * h) + lhs_.second(points[i]);
            matrix[i][i + 1] = 1 / (h * h) + lhs_.first(points[i]) / (2 * h);
            coefs[i] = rhs_(points[i]);
        }

        std::vector<double> solution;
        if constexpr (T == BoundaryConditionKind::FirstKind) {
            solution = matrix.SolveTridiagonalSystem(coefs);
        } else {
            solution = matrix.SolveSystem(coefs);
        }

        return solution;
    }

    std::vector<double> CalculateDiff(std::vector<double> const& first_solution,
                                      std::vector<double> const& second_solution) const {
        std::vector<double> diff(second_solution.size());
        for (std::size_t i = 0; i != first_solution.size(); ++i) {
            if (i % 2 == 0) {
                diff[i] = (second_solution[2 * i] - first_solution[i]) / 3;
            }
        }
        for (std::size_t i = 0; i != second_solution.size(); ++i) {
            if (i % 2 == 1) {
                diff[i] = (diff[i - 1] + diff[i + 1]) / 2;
            }
        }

        return diff;
    }

    double CalculateNorm(std::vector<double> const& first_solution,
                         std::vector<double> const& second_solution) const {
        std::vector<double> diff = CalculateDiff(first_solution, second_solution);
        double norm = 0;
        for (std::size_t i = 0; i != diff.size(); ++i) {
            norm = std::max(norm, std::abs(diff[i]));
        }
        return norm;
    }

    std::vector<double> SpecifySolution(std::vector<double> const& first_solution,
                                        std::vector<double> const& second_solution) const {
        std::vector<double> new_solution(second_solution.size());
        std::vector<double> diff = CalculateDiff(first_solution, second_solution);
        for (std::size_t i = 0; i != new_solution.size(); ++i) {
            new_solution[i] = second_solution[i] + diff[i];
        }

        return new_solution;
    }

    std::vector<double> CreateZeroSystemVector(double a, double b) const {
        model::Matrix m(2);
        if constexpr (T == BoundaryConditionKind::FirstKind) {
            m[0][0] = 1;
            m[0][1] = a;
            m[1][0] = 1;
            m[1][1] = b;
        } else if constexpr (T == BoundaryConditionKind::SecondKind) {
            m[0][0] = 0;
            m[0][1] = 1;
            m[1][0] = 0;
            m[1][1] = 1;
        } else if constexpr (T == BoundaryConditionKind::ThirdKind) {
            m[0][0] = boundary_condition_.left_value_coef;
            m[0][1] = boundary_condition_.left_value_coef * a +
                      boundary_condition_.left_derivative_coef;
            m[1][0] = boundary_condition_.right_value_coef;
            m[1][1] = boundary_condition_.right_value_coef * b +
                      boundary_condition_.right_derivative_coef;
        }
        if (m.Determinant() == 0) {
            throw std::runtime_error("Error: this task is not supported");
        }
        std::vector<double> vec = {boundary_condition_.left_value, boundary_condition_.right_value};
        std::vector<double> coefs0 = m.SolveSystem(vec);

        return coefs0;
    }

    std::pair<std::vector<double>, std::vector<double>> CreateFirstSystemVectors(double a,
                                                                                 double b) const {
        Matrix m1(2);
        std::vector<double> vec1(2);
        if constexpr (T == BoundaryConditionKind::FirstKind) {
            m1[0][0] = a;
            m1[0][1] = 1;
            m1[1][0] = b;
            m1[1][1] = 1;
            vec1[0] = -a * a;
            vec1[1] = -b * b;
        } else if constexpr (T == BoundaryConditionKind::SecondKind) {
            m1[0][0] = 1;
            m1[0][1] = 0;
            m1[1][0] = 1;
            m1[1][1] = 0;
            vec1[0] = -2 * a;
            vec1[1] = -2 * b;
        } else if constexpr (T == BoundaryConditionKind::ThirdKind) {
            m1[0][0] = a * boundary_condition_.left_value_coef +
                       boundary_condition_.left_derivative_coef;
            m1[0][1] = boundary_condition_.left_value_coef;
            m1[1][0] = b * boundary_condition_.right_value_coef +
                       boundary_condition_.right_derivative_coef;
            m1[1][1] = boundary_condition_.right_value_coef;
            vec1[0] = -boundary_condition_.left_value_coef * a * a -
                      2 * a * boundary_condition_.left_derivative_coef;
            vec1[1] = -boundary_condition_.right_value_coef * b * b -
                      2 * b * boundary_condition_.right_derivative_coef;
        }

        if (m1.Determinant() == 0) {
            throw std::runtime_error("Error: this task is not supported");
        }
        std::vector<double> coefs1 = m1.SolveSystem(vec1);

        std::vector<double> vec2(2);
        if constexpr (T == BoundaryConditionKind::FirstKind) {
            vec2[0] = -a * a * a;
            vec2[1] = -b * b * b;
        } else if constexpr (T == BoundaryConditionKind::SecondKind) {
            vec2[0] = -3 * a * a;
            vec2[1] = -3 * b * b;
        } else if constexpr (T == BoundaryConditionKind::ThirdKind) {
            vec2[0] = -boundary_condition_.left_value_coef * a * a * a -
                      3 * a * a * boundary_condition_.left_derivative_coef;
            vec2[1] = -boundary_condition_.right_value_coef * b * b * b -
                      3 * b * b * boundary_condition_.right_derivative_coef;
        }
        std::vector<double> coefs2 = m1.SolveSystem(vec2);

        return {std::move(coefs1), std::move(coefs2)};
    }

    std::vector<std::vector<double>> CalculateJacobiPolinomials(std::vector<double> const& points,
                                                                std::size_t n,
                                                                std::size_t k) const {
        if (n == 0) return {};
        std::vector<std::vector<double>> coefs(n, std::vector<double>(points.size()));
        for (std::size_t i = 0; i != points.size(); ++i) {
            coefs[0][i] = 1;
            if (n >= 2) {
                coefs[1][i] = (k + 1) * points[i];
            }
        }
        for (std::size_t i = 2; i < n; ++i) {
            for (std::size_t j = 0; j != points.size(); ++j) {
                coefs[i][j] =
                        ((i - 2 + k + 2) * (2 * (i - 2) + 2 * k + 3) * points[j] * coefs[i - 1][j] -
                         (i - 2 + k + 2) * (i - 2 + k + 1) * coefs[i - 2][j]) /
                        ((i - 2 + 2 * k + 2) * (i - 2 + 2));
            }
        }

        return coefs;
    }

    std::vector<std::vector<double>> CalculateJacobiDerivatives(std::vector<double> const& points,
                                                                std::size_t n,
                                                                std::size_t k) const {
        std::vector<std::vector<double>> coefs = CalculateJacobiPolinomials(points, n - 1, k + 1);
        std::vector<std::vector<double>> new_coefs(n, std::vector<double>(points.size()));
        for (std::size_t i = 1; i < n; ++i) {
            for (std::size_t j = 0; j != points.size(); ++j) {
                new_coefs[i][j] = (i + 2 * k + 1) * coefs[i - 1][j] / 2;
            }
        }

        return new_coefs;
    }

    std::vector<std::vector<double>> CalculateJacobiSecondDerivatives(
            std::vector<double> const& points, std::size_t n, std::size_t k) const {
        std::vector<std::vector<double>> coefs = CalculateJacobiDerivatives(points, n - 1, k + 1);
        std::vector<std::vector<double>> new_coefs(n, std::vector<double>(points.size()));
        for (std::size_t i = 1; i < n; ++i) {
            for (std::size_t j = 0; j != points.size(); ++j) {
                new_coefs[i][j] = (i + 2 * k + 1) * coefs[i - 1][j] / 2;
            }
        }

        return new_coefs;
    }

    double CalculateJacobiPolinomial(std::size_t n, std::size_t k, double x) const {
        if (n == 0) {
            return 1;
        }
        if (n == 1) {
            return (k + 1) * x;
        }

        double first_value = (n - 2 + k + 2) * (2 * (n - 2) + 2 * k + 3) * x *
                             CalculateJacobiPolinomial(n - 1, k, x);
        double second_value =
                (n - 2 + k + 2) * (n - 2 + k + 1) * CalculateJacobiPolinomial(n - 2, k, x);

        return (first_value - second_value) / ((n - 2 + 2 * k + 2) * (n - 2 + 2));
    }

    double CalculateJacobiDerivative(std::size_t n, std::size_t k, double x) const {
        if (n == 0) {
            return 0;
        }
        if (n == 1) {
            return k + 1;
        }

        return (n + 2 * k + 1) * CalculateJacobiPolinomial(n - 1, k + 1, x) / 2;
    }

    double CalculateJacobiSecondDerivative(std::size_t n, std::size_t k, double x) const {
        if (n == 0 || n == 1) {
            return 0;
        }

        return (n + 2 * k + 1) * CalculateJacobiDerivative(n - 1, k + 1, x) / 2;
    }

    Func CollocationMethod(std::vector<double> const& coefs0, std::vector<double> const& coefs1,
                           std::vector<double> const& coefs2, double a, double b,
                           std::size_t n) const {
        std::vector<double> points = CreateInnerPoints(a, b, n);
        std::vector<double> new_points(n);
        for (std::size_t i = 0; i != n; ++i) {
            new_points[i] = (2 * points[i] - b - a) / (b - a);
        }
        std::vector<std::vector<double>> values = CalculateJacobiPolinomials(new_points, n, 2);
        std::vector<std::vector<double>> first_derivatives =
                CalculateJacobiDerivatives(new_points, n, 2);
        std::vector<std::vector<double>> second_derivatives =
                CalculateJacobiSecondDerivatives(new_points, n, 2);

        auto jacobi_polinomial = [this](std::size_t i, double x) {
            return CalculateJacobiPolinomial(i, 2, x);
        };

        Matrix matrix(n);
        std::vector<double> vec(n);
        for (std::size_t i = 0; i != n; ++i) {
            for (std::size_t j = 0; j != n; ++j) {
                double second_derivative, first_derivative, value;
                if (j == 0) {
                    second_derivative = 2;
                    first_derivative = 2 * points[i] + coefs1[0];
                    value = points[i] * points[i] + coefs1[0] * points[i] + coefs1[1];
                } else if (j == 1) {
                    second_derivative = 6 * points[i];
                    first_derivative = 3 * points[i] * points[i] + coefs2[0];
                    value = points[i] * points[i] * points[i] + coefs2[0] * points[i] + coefs2[1];
                } else {
                    second_derivative =
                            (1 - new_points[i] * new_points[i]) *
                                    (1 - new_points[i] * new_points[i]) *
                                    second_derivatives[j - 2][i] * 4 / ((b - a) * (b - a)) -
                            first_derivatives[j - 2][i] * 8 * new_points[i] *
                                    (1 - new_points[i] * new_points[i]) * 4 / ((b - a) * (b - a)) +
                            values[j - 2][i] * 4 * (3 * new_points[i] * new_points[i] - 1) * 4 /
                                    ((b - a) * (b - a));
                    first_derivative = (1 - new_points[i] * new_points[i]) *
                                               (1 - new_points[i] * new_points[i]) *
                                               first_derivatives[j - 2][i] * 2 / (b - a) -
                                       values[j - 2][i] * 4 * new_points[i] *
                                               (1 - new_points[i] * new_points[i]) * 2 / (b - a);
                    value = (1 - new_points[i] * new_points[i]) *
                            (1 - new_points[i] * new_points[i]) * values[j - 2][i];
                }
                matrix[i][j] = second_derivative + lhs_.first(points[i]) * first_derivative +
                               lhs_.second(points[i]) * value;
            }
            vec[i] = rhs_(points[i]) - lhs_.first(points[i]) * coefs0[1] -
                     lhs_.second(points[i]) * (coefs0[0] + coefs0[1] * points[i]);
        }

        std::vector<double> solution = matrix.SolveSystem(vec);

        Func solution_func = [a, b, n, coefs0, coefs1, coefs2, solution,
                              jacobi_polinomial](double x) {
            if (x < a || x > b) return 0.0;
            double value = coefs0[0] + coefs0[1] * x;
            for (std::size_t i = 0; i != n; ++i) {
                if (i == 0) {
                    value += solution[i] * (x * x + coefs1[0] * x + coefs1[1]);
                } else if (i == 1) {
                    value += solution[i] * (x * x * x + coefs2[0] * x + coefs2[1]);
                } else {
                    double new_x = (2 * x - b - a) / (b - a);
                    value += solution[i] * (1 - new_x * new_x) * (1 - new_x * new_x) *
                             jacobi_polinomial(i - 2, new_x);
                }
            }
            return value;
        };

        return solution_func;
    }

    Func GalerkinMethod(std::vector<double> const& coefs0, std::vector<double> const& coefs1,
                        std::vector<double> const& coefs2, double a, double b,
                        std::size_t n) const {
        auto jacobi_polinomial = [this](std::size_t i, double x) {
            return CalculateJacobiPolinomial(i, 2, x);
        };

        Func f1 = [coefs1](double x) { return x * x + coefs1[0] * x + coefs1[1]; };
        Func f2 = [coefs2](double x) { return x * x * x + coefs2[0] * x + coefs2[1]; };
        std::vector<Func> jacobi_polinomials(n);
        for (std::size_t i = 0; i != n; ++i) {
            jacobi_polinomials[i] = [a, b, i, this](double x) {
                double new_x = (2 * x - b - a) / (b - a);
                return (1 - new_x * new_x) * (1 - new_x * new_x) *
                       CalculateJacobiPolinomial(i, 2, new_x);
            };
        }

        Func lf1 = [coefs1, this](double x) {
            return 2 + lhs_.first(x) * (2 * x + coefs1[0]) +
                   lhs_.second(x) * (x * x + coefs1[0] * x + coefs1[1]);
        };
        Func lf2 = [coefs2, this](double x) {
            return 6 * x + lhs_.first(x) * (3 * x * x + coefs2[0]) +
                   lhs_.second(x) * (x * x * x + coefs2[0] * x + coefs2[1]);
        };

        std::vector<Func> jacobi_derivatives(n);
        for (std::size_t i = 0; i != n; ++i) {
            jacobi_derivatives[i] = [a, b, i, this](double x) {
                double new_x = (2 * x - b - a) / (b - a);
                double jacobi_value = CalculateJacobiPolinomial(i, 2, new_x);
                double jacobi_derivative = CalculateJacobiDerivative(i, 2, new_x);
                double jacobi_second_derivative = CalculateJacobiSecondDerivative(i, 2, new_x);
                double second_derivative =
                        (1 - new_x * new_x) * (1 - new_x * new_x) * jacobi_second_derivative * 4 /
                                ((b - a) * (b - a)) -
                        jacobi_derivative * 8 * new_x * (1 - new_x * new_x) * 4 /
                                ((b - a) * (b - a)) +
                        jacobi_value * 4 * (3 * new_x * new_x - 1) * 4 / ((b - a) * (b - a));
                double first_derivative =
                        (1 - new_x * new_x) * (1 - new_x * new_x) * jacobi_derivative * 2 /
                                (b - a) -
                        jacobi_value * 4 * new_x * (1 - new_x * new_x) * 2 / (b - a);
                double value = (1 - new_x * new_x) * (1 - new_x * new_x) * jacobi_value;
                return second_derivative + lhs_.first(x) * first_derivative +
                       lhs_.second(x) * value;
            };
        }
        Func f0_derivative = [coefs0, this](double x) {
            return rhs_(x) - lhs_.first(x) * coefs0[1] -
                   lhs_.second(x) * (coefs0[0] + coefs0[1] * x);
        };

        Matrix matrix(n);
        std::vector<double> vec(n);
        for (std::size_t i = 0; i != n; ++i) {
            for (std::size_t j = 0; j != n; ++j) {
                Func coef_func;
                if (i == 0) {
                    if (j == 0) {
                        coef_func = [f1, lf1](double x) { return lf1(x) * f1(x); };
                    } else if (j == 1) {
                        coef_func = [f1, lf2](double x) { return lf2(x) * f1(x); };
                    } else {
                        coef_func = [f1, jacobi_derivatives, j](double x) {
                            return jacobi_derivatives[j - 2](x) * f1(x);
                        };
                    }
                } else if (i == 1) {
                    if (j == 0) {
                        coef_func = [f2, lf1](double x) { return lf1(x) * f2(x); };
                    } else if (j == 1) {
                        coef_func = [f2, lf2](double x) { return lf2(x) * f2(x); };
                    } else {
                        coef_func = [f2, jacobi_derivatives, j](double x) {
                            return jacobi_derivatives[j - 2](x) * f2(x);
                        };
                    }
                } else {
                    if (j == 0) {
                        coef_func = [jacobi_polinomials, lf1, i](double x) {
                            return lf1(x) * jacobi_polinomials[i - 2](x);
                        };
                    } else if (j == 1) {
                        coef_func = [jacobi_polinomials, lf2, i](double x) {
                            return lf2(x) * jacobi_polinomials[i - 2](x);
                        };
                    } else {
                        coef_func = [jacobi_polinomials, jacobi_derivatives, i, j](double x) {
                            return jacobi_derivatives[j - 2](x) * jacobi_polinomials[i - 2](x);
                        };
                    }
                }
                DefiniteIntegralCalculator calc{coef_func};
                matrix[i][j] = calc.Integral(a, b);
            }

            Func vec_func;
            if (i == 0) {
                vec_func = [f0_derivative, f1](double x) { return f0_derivative(x) * f1(x); };
            } else if (i == 1) {
                vec_func = [f0_derivative, f2](double x) { return f0_derivative(x) * f2(x); };
            } else {
                vec_func = [f0_derivative, jacobi_polinomials, i](double x) {
                    return f0_derivative(x) * jacobi_polinomials[i - 2](x);
                };
            }
            DefiniteIntegralCalculator vec_calc{vec_func};
            vec[i] = vec_calc.Integral(a, b);
        }

        std::vector<double> solution = matrix.SolveSystem(vec);

        Func solution_func = [a, b, n, coefs0, coefs1, coefs2, solution,
                              jacobi_polinomial](double x) {
            if (x < a || x > b) return 0.0;
            double value = coefs0[0] + coefs0[1] * x;
            for (std::size_t i = 0; i != n; ++i) {
                if (i == 0) {
                    value += solution[i] * (x * x + coefs1[0] * x + coefs1[1]);
                } else if (i == 1) {
                    value += solution[i] * (x * x * x + coefs2[0] * x + coefs2[1]);
                } else {
                    double new_x = (2 * x - b - a) / (b - a);
                    value += solution[i] * (1 - new_x * new_x) * (1 - new_x * new_x) *
                             jacobi_polinomial(i - 2, new_x);
                }
            }
            return value;
        };

        return solution_func;
    }

public:
    LinearODESolver() = default;

    LinearODESolver(std::pair<Func, Func> const& lhs, Func const& rhs,
                    BoundaryCondition<T> const& boundary_cond)
        : lhs_(lhs), rhs_(rhs) {
        if (boundary_cond.left_boundary >= boundary_cond.right_boundary) {
            throw std::invalid_argument("Error: the left boundary must be less than the right one");
        }
        if constexpr (T == BoundaryConditionKind::ThirdKind) {
            if (boundary_cond.left_value_coef == 0 && boundary_cond.left_derivative_coef == 0) {
                throw std::invalid_argument(
                        "Error: coefficients for the left boundary must not be both zeros");
            }
            if (boundary_cond.right_value_coef == 0 && boundary_cond.right_derivative_coef == 0) {
                throw std::invalid_argument(
                        "Error: coefficients for the right boundary must not be both zeros");
            }
            if (boundary_cond.left_value_coef * boundary_cond.left_derivative_coef > 0) {
                throw std::invalid_argument(
                        "Error: coefficients for the left boundary must have different signs");
            }
            if (boundary_cond.right_value_coef * boundary_cond.right_derivative_coef < 0) {
                throw std::invalid_argument(
                        "Error: coefficients for the right boundary must have the same sign");
            }
        }
        boundary_condition_ = boundary_cond;
    }

    LinearODESolver(std::pair<Func, Func>&& lhs, Func&& rhs, BoundaryCondition<T>&& boundary_cond)
        : lhs_(std::move(lhs)), rhs_(rhs) {
        if (boundary_cond.left_boundary >= boundary_cond.right_boundary) {
            throw std::invalid_argument("Error: the left boundary must be less than the right one");
        }
        if constexpr (T == BoundaryConditionKind::ThirdKind) {
            if (boundary_cond.left_value_coef == 0 && boundary_cond.left_derivative_coef == 0) {
                throw std::invalid_argument(
                        "Error: coefficients for the left boundary must not be both zeros");
            }
            if (boundary_cond.right_value_coef == 0 && boundary_cond.right_derivative_coef == 0) {
                throw std::invalid_argument(
                        "Error: coefficients for the right boundary must not be both zeros");
            }
            if (boundary_cond.left_value_coef * boundary_cond.left_derivative_coef > 0) {
                throw std::invalid_argument(
                        "Error: coefficients for the left boundary must have different signs");
            }
            if (boundary_cond.right_value_coef * boundary_cond.right_derivative_coef < 0) {
                throw std::invalid_argument(
                        "Error: coefficients for the right boundary must have the same sign");
            }
        }
        boundary_condition_ = std::move(boundary_cond);
    }

    std::pair<Func, SolveInfo> Solve(double eps) const {
        SolveInfo solve_info;
        double a = boundary_condition_.left_boundary;
        double b = boundary_condition_.right_boundary;
        std::size_t points_num = 5;
        double norm = 0;
        std::vector<double> first_solution, second_solution;
        first_solution = CalculateSolution(points_num);
        do {
            points_num *= 2;
            second_solution = CalculateSolution(points_num);
            norm = CalculateNorm(first_solution, second_solution);
            solve_info.emplace_back((b - a) / points_num, norm);
            if (norm > eps) {
                first_solution = std::move(second_solution);
            }
        } while (norm > eps);

        std::vector<double> final_solution = SpecifySolution(first_solution, second_solution);

        double h = (b - a) / points_num;
        std::vector<double> points = CreatePoints(a, b, h, points_num);
        Func solution_func = [a, b, final_solution, points](double x) {
            if (x < a || x > b) return 0.0;
            std::size_t index = std::lower_bound(points.begin(), points.end(), x) - points.begin();
            return final_solution[index];
        };

        return {std::move(solution_func), std::move(solve_info)};
    }

    Func Solve(std::size_t n, ODESolveMethod const method = ODESolveMethod::Collocation) const {
        double a = boundary_condition_.left_boundary;
        double b = boundary_condition_.right_boundary;
        std::vector<double> coefs0 = CreateZeroSystemVector(a, b);
        auto [coefs1, coefs2] = CreateFirstSystemVectors(a, b);

        Func solution;
        switch (method) {
            case ODESolveMethod::Collocation:
                solution = CollocationMethod(coefs0, coefs1, coefs2, a, b, n);
                break;
            case ODESolveMethod::Galerkin:
                solution = GalerkinMethod(coefs0, coefs1, coefs2, a, b, n);
                break;
            default:
                break;
        }

        return solution;
    }
};

}  // namespace model
