#include "task10.h"

#include <iostream>
#include <utility>

#include "optimization_methods.h"

namespace tasks {
void Semester6Task10() {
    std::cout << "Multivariable optimization" << '\n';
    std::cout << "Function: f(x) = (x - 2) * (x - 2)" << '\n';

    model::MultivariableFunc<1> f = [](std::array<double, 1> const& x) {
        return (x[0] - 2) * (x[0] - 2);
    };
    model::GradientFunc<1> f_gradient = [](std::array<double, 1> const& x) {
        std::array<double, 1> result = {2 * (x[0] - 2)};
        return result;
    };
    model::SecondDerivativeFunc<1> f_second_derivative = [](std::array<double, 1> const& x) {
        std::array<std::array<double, 1>, 1> result = {{2}};
        return result;
    };

    std::cout << "Precise local minimum: x0 = 2, f(x0) = 0" << '\n';
    std::cout << "Accuracy: 1e-5" << '\n';
    semester6_task10::CompareSolutions({2}, 0, f, f_gradient, f_second_derivative, 1e-5,
                                       {{0.1, 0.3}});

    std::cout << "Function: f(x, y) = ((x - 2) + (y - 1)) ^ 2 / 4 + ((x - 2) - (y - 1)) ^ 2"
              << '\n';

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
    model::SecondDerivativeFunc<2> f1_second_derivative = [](std::array<double, 2> const& x) {
        std::array<double, 2> first_row = {2.5, -1.5};
        std::array<double, 2> second_row = {-1.5, 2.5};
        std::array<std::array<double, 2>, 2> result = {std::move(first_row), std::move(second_row)};
        return result;
    };

    std::cout << "Precise local minimum: x0 = (2, 1), f(x0) = 0" << '\n';
    std::cout << "Accuracy: 1e-5" << '\n';
    semester6_task10::CompareSolutions(
            {2, 1}, 0, f1, f1_gradient, f1_second_derivative, 1e-5,
            {{0.1, 0.3}, {0.2, 0.3}, {0.3, 0.3}, {0.4, 0.3}, {0.5, 0.3}});

    std::cout << "Function: f(x, y) = 100 * (y - x^2)^2 + (1 - x)^2" << '\n';

    model::MultivariableFunc<2> f2 = [](std::array<double, 2> const& x) {
        double x1 = x[0];
        double x2 = x[1];
        return 100 * (x2 - x1 * x1) * (x2 - x1 * x1) + (1 - x1) * (1 - x1);
    };
    model::GradientFunc<2> f2_gradient = [](std::array<double, 2> const& x) {
        double x1 = x[0];
        double x2 = x[1];
        std::array<double, 2> result = {100 * 2 * (x2 - x1 * x1) * (-2 * x1) - 2 * (1 - x1),
                                        2 * (x2 - x1 * x1)};
        return result;
    };
    model::SecondDerivativeFunc<2> f2_second_derivative = [](std::array<double, 2> const& x) {
        double x1 = x[0];
        double x2 = x[1];
        std::array<double, 2> first_row = {-400 * (x2 - 3 * x1 * x1) + 2, -400 * x1};
        std::array<double, 2> second_row = {-4 * x1, 2};
        std::array<std::array<double, 2>, 2> result = {std::move(first_row), std::move(second_row)};
        return result;
    };

    std::cout << "Precise local minimum: x0 = (1, 1), f(x0) = 0" << '\n';
    std::cout << "Accuracy: 1e-5" << '\n';
    semester6_task10::CompareSolutions({1, 1}, 0, f2, f2_gradient, f2_second_derivative, 1e-5,
                                       {{0.1, 0.3}});
}

}  // namespace tasks
