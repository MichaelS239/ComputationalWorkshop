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

enum class BoundaryType { Equal, LessOrEqual };

template <std::size_t N>
struct Boundary {
    model::MultivariableFunc<N> boundary_func;
    std::array<model::MultivariableFunc<N>, N> boundary_derivatives;
    BoundaryType boundary_type;
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
MethodInfo<N> GradientProjectionMethod(model::MultivariableFunc<N> f,
                                       model::GradientFunc<N> gradient,
                                       std::vector<Boundary<N>> const& boundaries,
                                       std::array<double, N> const& x_start, double eps) {
    if (boundaries.size() >= N) {
        throw std::runtime_error("Error: too many boundaries");
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
                    coef += boundaries[j].boundary_derivatives[i](x0) *
                            boundaries[k].boundary_derivatives[i](x0);
                }
                matrix[k][j] = coef;
            }
            double vec_coef = 0;
            for (std::size_t i = 0; i != N; ++i) {
                vec_coef -= boundaries[k].boundary_derivatives[i](x0) * gradient_value[i];
            }
            vec[k] = vec_coef;
        }
        std::vector<double> lambdas = matrix.SolveSystem(vec);
        std::vector<double> gradient_projection(N);
        for (std::size_t i = 0; i != N; ++i) {
            gradient_projection[i] = gradient_value[i];
            for (std::size_t j = 0; j != boundaries.size(); ++j) {
                gradient_projection[i] += lambdas[j] * boundaries[j].boundary_derivatives[i](x0);
            }
        }
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
        } while (dif >= 0);
        x0 = x1;
    }
    double min_value = f(x0);

    return {std::move(x0), min_value, iter_num};
}

template <std::size_t N>
MethodInfo<N> HeavyBallMethod(model::MultivariableFunc<N> f, model::GradientFunc<N> gradient,
                              std::array<double, N> const& start_x, double eps, double alpha = 0.1,
                              double beta = 0.3) {
    std::array<double, N> x0;
    std::array<double, N> x1;
    std::array<double, N> x2 = start_x;
    std::array<double, N> gradient_value;
    for (std::size_t i = 0; i != N; ++i) {
        x0[i] = 0;
        x1[i] = 0;
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
MethodInfo<N> PenaltyFunctionMethod(model::MultivariableFunc<N> f, model::GradientFunc<N> gradient,
                                    std::vector<Boundary<N>> const& boundaries,
                                    std::array<double, N> const& x_start, double eps) {
    if (boundaries.size() >= N) {
        throw std::runtime_error("Error: too many boundaries");
    }
    model::MultivariableFunc<N> penalty_func = [boundaries](std::array<double, N> const& x) {
        double sum = 0;
        for (std::size_t i = 0; i != boundaries.size(); ++i) {
            if (boundaries[i].boundary_type == BoundaryType::LessOrEqual) {
                double max_value = std::max(0.0, boundaries[i].boundary_func(x));
                sum += max_value * max_value;
            } else if (boundaries[i].boundary_type == BoundaryType::Equal) {
                double max_value = std::abs(boundaries[i].boundary_func(x));
                sum += max_value * max_value;
            }
        }
        return sum;
    };

    model::GradientFunc<N> penalty_func_gradient = [boundaries](std::array<double, N> const& x) {
        std::array<double, N> gradient;
        for (std::size_t i = 0; i != N; ++i) {
            gradient[i] = 0;
        }
        for (std::size_t i = 0; i != boundaries.size(); ++i) {
            if (boundaries[i].boundary_type == BoundaryType::LessOrEqual) {
                if (boundaries[i].boundary_func(x) > 0) {
                    for (std::size_t j = 0; j != N; ++j) {
                        gradient[j] += 2 * boundaries[i].boundary_func(x) *
                                       boundaries[i].boundary_derivatives[j](x);
                    }
                }
            } else if (boundaries[i].boundary_type == BoundaryType::Equal) {
                for (std::size_t j = 0; j != N; ++j) {
                    gradient[j] += 2 * boundaries[i].boundary_func(x) *
                                   boundaries[i].boundary_derivatives[j](x);
                }
            }
        }
        return gradient;
    };

    double diff = 2 * eps;
    double alpha = 2;
    std::array<double, N> x0 = x_start;
    std::size_t iter_num = 0;
    while (diff > eps) {
        if (iter_num == 200) {
            break;
        }
        ++iter_num;
        alpha *= alpha;
        model::MultivariableFunc<N> teta_func = [f, penalty_func,
                                                 alpha](std::array<double, N> const& x) {
            return f(x) + alpha * penalty_func(x);
        };

        model::GradientFunc<N> teta_gradient = [gradient, penalty_func_gradient,
                                                alpha](std::array<double, N> const& x) {
            std::array<double, N> result = gradient(x);
            std::array<double, N> penalty = penalty_func_gradient(x);
            for (std::size_t i = 0; i != N; ++i) {
                result[i] += alpha * penalty[i];
            }
            return result;
        };

        MethodInfo<N> minimum = HeavyBallMethod(teta_func, teta_gradient, x0, eps);
        x0 = minimum.local_minimum;
        diff = alpha * penalty_func(x0);
    }
    double min_value = f(x0);

    return {std::move(x0), min_value, iter_num};
}

template <std::size_t N>
MethodInfo<N> BarrierFunctionMethod(model::MultivariableFunc<N> f, model::GradientFunc<N> gradient,
                                    std::vector<Boundary<N>> const& boundaries,
                                    std::array<double, N> const& x_start, double eps) {
    if (boundaries.size() >= N) {
        throw std::runtime_error("Error: too many boundaries");
    }
    model::MultivariableFunc<N> barrier_func = [boundaries](std::array<double, N> const& x) {
        double sum = 0;
        for (std::size_t i = 0; i != boundaries.size(); ++i) {
            sum -= 1 / boundaries[i].boundary_func(x);
        }
        return sum;
    };

    model::GradientFunc<N> barrier_func_gradient = [boundaries](std::array<double, N> const& x) {
        std::array<double, N> gradient;
        for (std::size_t i = 0; i != N; ++i) {
            gradient[i] = 0;
        }
        for (std::size_t i = 0; i != boundaries.size(); ++i) {
            for (std::size_t j = 0; j != N; ++j) {
                gradient[j] += 1 /
                               (boundaries[i].boundary_func(x) * boundaries[i].boundary_func(x)) *
                               boundaries[i].boundary_derivatives[j](x);
            }
        }
        return gradient;
    };

    double diff = 2 * eps;
    double mu = 0.5;
    std::array<double, N> x0 = x_start;
    std::size_t iter_num = 0;
    while (diff > eps) {
        if (iter_num == 200) {
            break;
        }
        ++iter_num;
        mu *= mu;
        model::MultivariableFunc<N> teta_func = [f, barrier_func,
                                                 mu](std::array<double, N> const& x) {
            return f(x) + mu * barrier_func(x);
        };

        model::GradientFunc<N> teta_gradient = [gradient, barrier_func_gradient,
                                                mu](std::array<double, N> const& x) {
            std::array<double, N> result = gradient(x);
            std::array<double, N> barrier = barrier_func_gradient(x);
            for (std::size_t i = 0; i != N; ++i) {
                result[i] += mu * barrier[i];
            }
            return result;
        };

        MethodInfo<N> minimum = HeavyBallMethod(teta_func, teta_gradient, x0, eps);
        x0 = minimum.local_minimum;
        diff = mu * barrier_func(x0);
    }
    double min_value = f(x0);

    return {std::move(x0), min_value, iter_num};
}

}  // namespace semester6_task11
