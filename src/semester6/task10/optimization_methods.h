#pragma once

#include <array>
#include <cmath>
#include <cstddef>

#include "model/function.h"

namespace semester6_task10 {

template <std::size_t N>
struct MethodInfo {
    std::array<double, N> local_minimum;
    double min_value;
    std::size_t iteration_number;
};

template <std::size_t N>
double Norm(std::array<double, N> const& vector) {
    double norm = 0;
    for (std::size_t i = 0; i != N; ++i) {
        norm += vector[i] * vector[i];
    }
    return std::sqrt(norm);
}

template <std::size_t N>
MethodInfo<N> GradientDescent(model::MultivariableFunc<N> f, model::GradientFunc<N> gradient,
                              double eps, double gamma = 0.1) {
    std::array<double, N> x0;
    std::array<double, N> gradient_value;
    for (std::size_t i = 0; i != N; ++i) {
        x0[i] = 0;
        gradient_value[i] = 0;
    }
    double norm = 0;
    std::size_t iter_num = 0;
    do {
        if (iter_num == 200) {
            break;
        }
        ++iter_num;
        for (std::size_t i = 0; i != N; ++i) {
            x0[i] -= gamma * gradient_value[i];
        }
        gradient_value = gradient(x0);
        norm = Norm(gradient_value);
    } while (norm >= eps);
    double min_value = f(x0);

    return {std::move(x0), min_value, iter_num};
}

template <std::size_t N>
MethodInfo<N> HeavyBall(model::MultivariableFunc<N> f, model::GradientFunc<N> gradient, double eps,
                        double alpha = 0.1, double beta = 0.1) {
    std::array<double, N> x0;
    std::array<double, N> x1;
    std::array<double, N> x2;
    std::array<double, N> gradient_value;
    for (std::size_t i = 0; i != N; ++i) {
        x0[i] = 0;
        x1[i] = 0;
        x2[i] = 0;
        gradient_value[i] = 0;
    }
    double norm = 0;
    std::size_t iter_num = 0;
    do {
        if (iter_num == 200) {
            break;
        }
        ++iter_num;
        for (std::size_t i = 0; i != N; ++i) {
            x2[i] = x2[i] - alpha * gradient_value[i] + beta * (x1[i] - x0[i]);
        }
        gradient_value = gradient(x2);
        norm = Norm(gradient_value);
        x0 = x1;
        x1 = x2;
    } while (norm >= eps);
    double min_value = f(x2);

    return {std::move(x2), min_value, iter_num};
}

}  // namespace semester6_task10
