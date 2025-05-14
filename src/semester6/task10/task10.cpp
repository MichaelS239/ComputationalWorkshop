#include "task10.h"

#include <iostream>

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

    std::cout << "Precise local minimum: x0 = 2, f(x0) = 0" << '\n';
    std::cout << "Accuracy: 1e-5" << '\n';
    std::cout << "Speed: 0.1" << '\n';

    semester6_task10::MethodInfo gradient_info =
            semester6_task10::GradientDescent(f, f_gradient, 1e-5, 0.1);
    std::cout << "Gradient descent: x0 = " << gradient_info.local_minimum[0]
              << ", f(x0) = " << gradient_info.min_value
              << ", iteration number: " << gradient_info.iteration_number << '\n';
    semester6_task10::MethodInfo heavy_ball_info =
            semester6_task10::HeavyBall(f, f_gradient, 1e-5, 0.1);
    std::cout << "Heavy ball: x0 = " << heavy_ball_info.local_minimum[0]
              << ", f(x0) = " << heavy_ball_info.min_value
              << ", iteration number: " << heavy_ball_info.iteration_number << '\n';
    semester6_task10::CompareSolutions({2}, 0, gradient_info, heavy_ball_info);

    std::cout << "Function: f(x) = ((x - 2) + (y - 1)) ^ 2 / 4 + ((x - 2) - (y - 1)) ^ 2" << '\n';

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

    std::cout << "Precise local minimum: x0 = (2, 1), f(x0) = 0" << '\n';
    std::cout << "Accuracy: 1e-5" << '\n';

    std::cout << "Speed: 0.1" << '\n';
    semester6_task10::MethodInfo gradient_info1 =
            semester6_task10::GradientDescent(f1, f1_gradient, 1e-5, 0.1);
    std::cout << "Gradient descent: f(x0) = " << gradient_info1.min_value
              << ", iteration number: " << gradient_info1.iteration_number << '\n';
    semester6_task10::MethodInfo heavy_ball_info1 =
            semester6_task10::HeavyBall(f1, f1_gradient, 1e-5, 0.1);
    std::cout << "Heavy ball: f(x0) = " << heavy_ball_info1.min_value
              << ", iteration number: " << heavy_ball_info1.iteration_number << '\n';
    semester6_task10::CompareSolutions({2, 1}, 0, gradient_info1, heavy_ball_info1);

    std::cout << "Speed: 0.45" << '\n';
    semester6_task10::MethodInfo gradient_info2 =
            semester6_task10::GradientDescent(f1, f1_gradient, 1e-5, 0.45);
    std::cout << "Gradient descent: f(x0) = " << gradient_info2.min_value
              << ", iteration number: " << gradient_info2.iteration_number << '\n';
    semester6_task10::MethodInfo heavy_ball_info2 =
            semester6_task10::HeavyBall(f1, f1_gradient, 1e-5, 0.45);
    std::cout << "Heavy ball: f(x0) = " << heavy_ball_info2.min_value
              << ", iteration number: " << heavy_ball_info2.iteration_number << '\n';
    semester6_task10::CompareSolutions({2, 1}, 0, gradient_info2, heavy_ball_info2);

    std::cout << "Speed: 0.5" << '\n';
    semester6_task10::MethodInfo gradient_info3 =
            semester6_task10::GradientDescent(f1, f1_gradient, 1e-5, 0.5);
    std::cout << "Gradient descent: f(x0) = " << gradient_info3.min_value
              << ", iteration number: " << gradient_info3.iteration_number << '\n';
    semester6_task10::MethodInfo heavy_ball_info3 =
            semester6_task10::HeavyBall(f1, f1_gradient, 1e-5, 0.5);
    std::cout << "Heavy ball: f(x0) = " << heavy_ball_info3.min_value
              << ", iteration number: " << heavy_ball_info3.iteration_number << '\n';
    semester6_task10::CompareSolutions({2, 1}, 0, gradient_info3, heavy_ball_info3);

    /*model::MultivariableFunc<2> f = [](std::array<double, 2> const& x) {
        double x1 = x[0];
        double x2 = x[1];
        return 100 * (x2 - x1 * x1) * (x2 - x1 * x1) + (1 - x1) * (1 - x1);
    };*/
}

}  // namespace tasks
