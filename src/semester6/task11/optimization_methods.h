#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <vector_util.h>

#include "model/function.h"
#include "model/matrix.h"

namespace semester6_task11 {

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
MethodInfo<N> HeavyBallMethod(model::MultivariableFunc<N> f, model::GradientFunc<N> gradient,
                              double eps, double alpha = 0.1, double beta = 0.3) {
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

template <std::size_t N>
MethodInfo<N> GradientProjectionMethod(
        model::MultivariableFunc<N> f, model::GradientFunc<N> gradient,
        std::vector<model::MultivariableFunc<N>> const& boundaries,
        std::vector<std::vector<model::MultivariableFunc<N>>> const& boundaries_derivatives,
        std::array<double, N> const& x_start, double eps) {
    if (boundaries.size() > N) {
        throw std::runtime_error("Error: too many boundaries");
    }
    if (boundaries_derivatives.size() != boundaries.size()) {
        throw std::runtime_error(
                "Error: the number of boundaries should be equal to the number of their "
                "derivatives");
    }
    double h = 0.1;
    std::array<double, N> x0 = x_start;
    std::array<double, N> gradient_value;
    for (std::size_t i = 0; i != N; ++i) {
        gradient_value[i] = 0;
    }
    double norm = 2 * eps;
    std::size_t iter_num = 0;
    while (norm >= eps) {
        if (iter_num == 200) {
            break;
        }
        ++iter_num;
        gradient_value = gradient(x0);
        model::Matrix matrix(boundaries.size());
        std::vector<double> vec(boundaries.size());
        for (std::size_t k = 0; k != boundaries.size(); ++k) {
            for (std::size_t j = 0; j != boundaries.size(); ++j) {
                double coef = 0;
                for (std::size_t i = 0; i != N; ++i) {
                    coef += boundaries_derivatives[j][i](x0) * boundaries_derivatives[k][i](x0);
                }
                matrix[k][j] = coef;
            }
            double vec_coef = 0;
            for (std::size_t i = 0; i != boundaries.size(); ++i) {
                vec_coef -= boundaries_derivatives[k][i](x0) * gradient_value[i];
            }
            vec[k] = vec_coef;
        }
        // std::cout << matrix.ToString() << '\n';
        // std::cout << vec[0] << '\n';
        std::vector<double> lambdas = matrix.SolveSystem(vec);
        /*for (std::size_t i = 0; i != boundaries.size(); ++i) {
            std::cout << lambdas[i] << ' ';
        }
        std::cout << '\n';*/
        std::vector<double> gradient_projection(N);
        for (std::size_t i = 0; i != N; ++i) {
            gradient_projection[i] = gradient_value[i];
            for (std::size_t j = 0; j != boundaries.size(); ++j) {
                gradient_projection[i] += lambdas[j] * boundaries_derivatives[j][i](x0);
            }
        }
        /*for (std::size_t i = 0; i != N; ++i) {
            std::cout << gradient_projection[i] << ' ';
        }
        std::cout << '\n';*/
        norm = util::Norm(gradient_projection);
        std::vector<double> directions(N);
        for (std::size_t i = 0; i != N; ++i) {
            directions[i] = gradient_projection[i] / norm;
        }
        std::array<double, N> x1;
        double dif = 0;
        h *= 2;
        do {
            h /= 2;
            for (std::size_t i = 0; i != N; ++i) {
                x1[i] = x0[i] - h * directions[i];
            }
            dif = f(x1) - f(x0);
            /*for (std::size_t i = 0; i != N; ++i) {
                std::cout << x1[i] << ' ';
            }
            std::cout << '\n';*/
        } while (dif >= 0);
        x0 = x1;
        /*for (std::size_t i = 0; i != N; ++i) {
            std::cout << x0[i] << ' ';
        }
        std::cout << '\n';*/
    }
    double min_value = f(x0);

    return {std::move(x0), min_value, iter_num};
}

}  // namespace semester6_task11
