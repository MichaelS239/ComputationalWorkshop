#include "task11.h"

#include <iostream>
#include <utility>

#include "optimization_methods.h"

namespace tasks {
void Semester6Task11() {
    std::cout << "Multivariable optimization with constraints" << '\n';
    std::cout << '\n';
    std::cout << "Function: f(x, y) = ((x - 2) + (y - 1)) ^ 2 / 4 + ((x - 2) - (y - 1)) ^ 2"
              << '\n';
    std::cout << "Precise local minimum: x0 = (2, 1), f(x0) = 0" << '\n';
    std::cout << "Accuracy: 1e-5" << '\n';

    model::MultivariableFunc<2> f1 = [](std::array<double, 2> const& x) {
        return ((x[0] - 2) + (x[1] - 1)) * ((x[0] - 2) + (x[1] - 1)) / 4 +
               ((x[0] - 2) - (x[1] - 1)) * ((x[0] - 2) - (x[1] - 1));
    };
    model::GradientFunc<2> f1_gradient = [](std::array<double, 2> const& x) {
        std::array<double, 2> result = {
                2 * ((x[0] - 2) + (x[1] - 1)) / 4 + 2 * ((x[0] - 2) - (x[1] - 1)),
                2 * ((x[0] - 2) + (x[1] - 1)) / 4 - 2 * ((x[0] - 2) - (x[1] - 1))};
        return result;
    };
    std::cout << "D = {(x, y) : x = 2}" << '\n';
    model::MultivariableFunc<2> f1_boundary = [](std::array<double, 2> const& x) {
        return x[0] - 2;
    };
    model::MultivariableFunc<2> f1_boundary_derivative_x = [](std::array<double, 2> const& x) {
        return 1;
    };
    model::MultivariableFunc<2> f1_boundary_derivative_y = [](std::array<double, 2> const& x) {
        return 0;
    };

    semester6_task11::CompareProjectionAndPenaltySolutions(
            {2, 1}, 0, f1, f1_gradient,
            {{f1_boundary,
              {f1_boundary_derivative_x, f1_boundary_derivative_y},
              semester6_task11::BoundaryType::Equal}},
            {2, 0}, {0.5, 0.5}, 1e-5);
    std::cout << "D = {(x, y) : y = 1}" << '\n';
    model::MultivariableFunc<2> f1_boundary1 = [](std::array<double, 2> const& x) {
        return x[1] - 1;
    };
    model::MultivariableFunc<2> f1_boundary1_derivative_x = [](std::array<double, 2> const& x) {
        return 0;
    };
    model::MultivariableFunc<2> f1_boundary1_derivative_y = [](std::array<double, 2> const& x) {
        return 1;
    };
    semester6_task11::CompareProjectionAndPenaltySolutions(
            {2, 1}, 0, f1, f1_gradient,
            {{f1_boundary1,
              {f1_boundary1_derivative_x, f1_boundary1_derivative_y},
              semester6_task11::BoundaryType::Equal}},
            {0, 1}, {0.5, 0.5}, 1e-5);

    std::cout << "D = {(x, y) : y = 0.5 * x}" << '\n';
    model::MultivariableFunc<2> f1_boundary4 = [](std::array<double, 2> const& x) {
        return x[1] - 0.5 * x[0];
    };
    model::MultivariableFunc<2> f1_boundary4_derivative_x = [](std::array<double, 2> const& x) {
        return -0.5;
    };
    model::MultivariableFunc<2> f1_boundary4_derivative_y = [](std::array<double, 2> const& x) {
        return 1;
    };
    semester6_task11::CompareProjectionAndPenaltySolutions(
            {2, 1}, 0, f1, f1_gradient,
            {{f1_boundary4,
              {f1_boundary4_derivative_x, f1_boundary4_derivative_y},
              semester6_task11::BoundaryType::Equal}},
            {0, 0}, {0.5, 0.5}, 1e-5);

    std::cout << "D = {(x, y) : y <= 0.5 * x}" << '\n';
    model::MultivariableFunc<2> f1_boundary5 = [](std::array<double, 2> const& x) {
        return x[1] - 0.5 * x[0];
    };
    model::MultivariableFunc<2> f1_boundary5_derivative_x = [](std::array<double, 2> const& x) {
        return -0.5;
    };
    model::MultivariableFunc<2> f1_boundary5_derivative_y = [](std::array<double, 2> const& x) {
        return 1;
    };
    semester6_task11::ComparePenaltyAndBarrierSolutions(
            {2, 1}, 0, f1, f1_gradient,
            {{f1_boundary5,
              {f1_boundary5_derivative_x, f1_boundary5_derivative_y},
              semester6_task11::BoundaryType::LessOrEqual}},
            {0.5, 0.5}, {0.5, 0}, 1e-5);

    std::cout << "D = {(x, y) : y <= x}" << '\n';
    model::MultivariableFunc<2> f1_boundary6 = [](std::array<double, 2> const& x) {
        return x[1] - x[0];
    };
    model::MultivariableFunc<2> f1_boundary6_derivative_x = [](std::array<double, 2> const& x) {
        return -1;
    };
    model::MultivariableFunc<2> f1_boundary6_derivative_y = [](std::array<double, 2> const& x) {
        return 1;
    };
    semester6_task11::ComparePenaltyAndBarrierSolutions(
            {2, 1}, 0, f1, f1_gradient,
            {{f1_boundary6,
              {f1_boundary6_derivative_x, f1_boundary6_derivative_y},
              semester6_task11::BoundaryType::LessOrEqual}},
            {0, 0.5}, {0.5, 0}, 1e-5);

    std::cout << "D = {(x, y) : y <= 2}" << '\n';
    model::MultivariableFunc<2> f1_boundary7 = [](std::array<double, 2> const& x) {
        return x[1] - 2;
    };
    model::MultivariableFunc<2> f1_boundary7_derivative_x = [](std::array<double, 2> const& x) {
        return 0;
    };
    model::MultivariableFunc<2> f1_boundary7_derivative_y = [](std::array<double, 2> const& x) {
        return 1;
    };
    semester6_task11::ComparePenaltyAndBarrierSolutions(
            {2, 1}, 0, f1, f1_gradient,
            {{f1_boundary7,
              {f1_boundary7_derivative_x, f1_boundary7_derivative_y},
              semester6_task11::BoundaryType::LessOrEqual}},
            {0, 3}, {0, 0}, 1e-5);

    std::cout << "D = {(x, y) : x^2 + y^2 <= 5}, x0 = (-1, -2)" << '\n';
    model::MultivariableFunc<2> f1_boundary8 = [](std::array<double, 2> const& x) {
        return x[0] * x[0] + x[1] * x[1] - 5;
    };
    model::MultivariableFunc<2> f1_boundary8_derivative_x = [](std::array<double, 2> const& x) {
        return 2 * x[0];
    };
    model::MultivariableFunc<2> f1_boundary8_derivative_y = [](std::array<double, 2> const& x) {
        return 2 * x[1];
    };
    semester6_task11::ComparePenaltyAndBarrierSolutions(
            {2, 1}, 0, f1, f1_gradient,
            {{f1_boundary8,
              {f1_boundary8_derivative_x, f1_boundary8_derivative_y},
              semester6_task11::BoundaryType::LessOrEqual}},
            {2, 2}, {0, 0}, 1e-5);
}

}  // namespace tasks
