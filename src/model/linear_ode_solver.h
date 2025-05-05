#pragma once

#include <algorithm>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "boundary_condition.h"
#include "matrix.h"

namespace model {

using SolveInfo = std::vector<std::pair<double, double>>;

template <BoundaryConditionKind T>
class LinearODESolver {
private:
    std::pair<std::function<double(double)>, std::function<double(double)>> lhs_;
    std::function<double(double)> rhs_;
    BoundaryCondition<T> boundary_condition_;

    std::vector<double> CreatePoints(double a, double b, double h, std::size_t n) const {
        std::vector<double> points(n + 1);
        for (std::size_t i = 0; i != n + 1; ++i) {
            points[i] = a + i * h;
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

public:
    LinearODESolver() = default;

    LinearODESolver(
            std::pair<std::function<double(double)>, std::function<double(double)>> const& lhs,
            std::function<double(double)> const& rhs, BoundaryCondition<T> const& boundary_cond)
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

    LinearODESolver(std::pair<std::function<double(double)>, std::function<double(double)>>&& lhs,
                    std::function<double(double)>&& rhs, BoundaryCondition<T>&& boundary_cond)
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

    std::pair<std::function<double(double)>, SolveInfo> Solve(double eps) const {
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
        std::function<double(double)> solution_func = [a, b, final_solution, points](double x) {
            if (x < a || x > b) return 0.0;
            std::size_t index = std::lower_bound(points.begin(), points.end(), x) - points.begin();
            return final_solution[index];
        };

        return {std::move(solution_func), std::move(solve_info)};
    }
};

}  // namespace model
