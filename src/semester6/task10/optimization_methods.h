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
                              double eps) {
    std::array<double, N> x0;
    std::array<double, N> gradient_value;
    double norm;
    double gamma = 0.1;
    std::size_t iter_num = 0;
    do {
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

}  // namespace semester6_task10
