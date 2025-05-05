#include "linear_ode_solver.h"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <utility>

#include "matrix.h"

namespace model {

std::vector<double> LinearODESolver::CreatePoints(double a, double b, double h,
                                                  std::size_t n) const {
    std::vector<double> points(n + 1);
    for (std::size_t i = 0; i != n + 1; ++i) {
        points[i] = a + i * h;
    }
    return points;
}

std::vector<double> LinearODESolver::CalculateSolution(std::size_t points_num) const {
    double a = boundary_condition_.left_boundary;
    double b = boundary_condition_.right_boundary;
    double h = (b - a) / points_num;
    std::vector<double> points = CreatePoints(a, b, h, points_num);
    Matrix matrix(points_num + 1);
    std::vector<double> coefs(points_num + 1);
    if (boundary_condition_.kind == BoundaryConditionKind::FirstKind) {
        matrix[0][0] = 1;
        coefs[0] = boundary_condition_.left_value;
        matrix[points_num][points_num] = 1;
        coefs[points_num] = boundary_condition_.right_value;
    } else if (boundary_condition_.kind == BoundaryConditionKind::SecondKind) {
        matrix[0][0] = -3 / (2 * h);
        matrix[0][1] = 2 / h;
        matrix[0][2] = -1 / (2 * h);
        coefs[0] = boundary_condition_.left_value;
        matrix[points_num][points_num] = 3 / (2 * h);
        matrix[points_num][points_num - 1] = -2 / h;
        matrix[points_num][points_num - 2] = 1 / (2 * h);
        coefs[points_num] = boundary_condition_.right_value;
    }
    for (std::size_t i = 1; i != points_num; ++i) {
        matrix[i][i - 1] = 1 / (h * h) - lhs_.first(points[i]) / (2 * h);
        matrix[i][i] = -2 / (h * h) + lhs_.second(points[i]);
        matrix[i][i + 1] = 1 / (h * h) + lhs_.first(points[i]) / (2 * h);
        coefs[i] = rhs_(points[i]);
    }

    std::vector<double> solution;
    if (boundary_condition_.kind == BoundaryConditionKind::FirstKind) {
        solution = matrix.SolveTridiagonalSystem(coefs);
    } else {
        solution = matrix.SolveSystem(coefs);
    }

    return solution;
}

std::vector<double> LinearODESolver::CalculateDiff(
        std::vector<double> const& first_solution,
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

double LinearODESolver::CalculateNorm(std::vector<double> const& first_solution,
                                      std::vector<double> const& second_solution) const {
    std::vector<double> diff = CalculateDiff(first_solution, second_solution);
    double norm = 0;
    for (std::size_t i = 0; i != diff.size(); ++i) {
        norm = std::max(norm, std::abs(diff[i]));
    }
    return norm;
}

std::vector<double> LinearODESolver::SpecifySolution(
        std::vector<double> const& first_solution,
        std::vector<double> const& second_solution) const {
    std::vector<double> new_solution(second_solution.size());
    std::vector<double> diff = CalculateDiff(first_solution, second_solution);
    for (std::size_t i = 0; i != new_solution.size(); ++i) {
        new_solution[i] = second_solution[i] + diff[i];
    }

    return new_solution;
}

std::pair<std::function<double(double)>, SolveInfo> LinearODESolver::Solve(double eps) const {
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

}  // namespace model
