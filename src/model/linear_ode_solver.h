#pragma once

#include <cstddef>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "boundary_condition.h"

namespace model {

using SolveInfo = std::vector<std::pair<double, double>>;

class LinearODESolver {
private:
    std::pair<std::function<double(double)>, std::function<double(double)>> lhs_;
    std::function<double(double)> rhs_;
    BoundaryCondition boundary_condition_;

    std::vector<double> CreatePoints(double a, double b, double h, std::size_t n) const;
    std::vector<double> CalculateSolution(std::size_t points_num) const;
    std::vector<double> CalculateDiff(std::vector<double> const& first_solution,
                                      std::vector<double> const& second_solution) const;
    double CalculateNorm(std::vector<double> const& first_solution,
                         std::vector<double> const& second_solution) const;
    std::vector<double> SpecifySolution(std::vector<double> const& first_solution,
                                        std::vector<double> const& second_solution) const;

public:
    LinearODESolver() = default;

    LinearODESolver(
            std::pair<std::function<double(double)>, std::function<double(double)>> const& lhs,
            std::function<double(double)> const& rhs, BoundaryCondition const& boundary_cond)
        : lhs_(lhs), rhs_(rhs) {
        if (boundary_cond.left_boundary >= boundary_cond.right_boundary) {
            throw std::invalid_argument("Error: the left boundary must be less than the right one");
        }
        boundary_condition_ = boundary_cond;
    }

    LinearODESolver(std::pair<std::function<double(double)>, std::function<double(double)>>&& lhs,
                    std::function<double(double)>&& rhs, BoundaryCondition&& boundary_cond)
        : lhs_(std::move(lhs)), rhs_(rhs) {
        if (boundary_cond.left_boundary >= boundary_cond.right_boundary) {
            throw std::invalid_argument("Error: the left boundary must be less than the right one");
        }
        boundary_condition_ = std::move(boundary_cond);
    }

    std::pair<std::function<double(double)>, SolveInfo> Solve(double eps) const;
};

}  // namespace model
