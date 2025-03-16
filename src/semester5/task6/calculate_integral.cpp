#include "calculate_integral.h"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

#include "model/polynomial.h"

namespace semester5_task6 {

double f(double x) {
    return std::sin(x);
}

std::pair<std::vector<std::vector<double>>, std::vector<double>> CreatePointsSystem(
        std::size_t num_points, std::vector<double> const& moments) {
    std::vector<std::vector<double>> matrix =
            std::vector<std::vector<double>>(num_points, std::vector<double>(num_points));
    std::vector<double> vector = std::vector<double>(num_points);
    for (std::size_t i = 0; i != matrix.size(); ++i) {
        for (std::size_t j = 0; j != matrix.size(); ++j) {
            matrix[i][j] = moments[i + j];
        }
    }

    for (std::size_t i = 0; i != vector.size(); ++i) {
        vector[i] = -moments[i + num_points];
    }

    return {matrix, vector};
}

std::vector<double> CalculateMoments(std::size_t n) {
    std::vector<double> moments(n);
    for (std::size_t i = 0; i != n; ++i) {
        moments[i] = 1.0 / ((i + 2) * (i + 2));
    }
    return moments;
}

std::vector<double> CalculateGaussMoments(std::size_t n) {
    std::vector<double> moments(n);
    for (std::size_t i = 0; i != n; ++i) {
        if (i % 2 == 0) {
            moments[i] = (double)2 / (i + 1);
        }
    }
    return moments;
}

std::vector<double> CalculateMellerMoments(std::size_t n) {
    std::vector<double> moments(n);
    for (std::size_t i = 0; i != n; ++i) {
        if (i == 0) {
            moments[i] = std::numbers::pi;
        } else if (i % 2 == 0) {
            moments[i] = 2 * std::tgamma(0.5) * std::tgamma(i / 2 + 0.5) / (i * std::tgamma(i / 2));
        }
    }
    return moments;
}

}  // namespace semester5_task6
