#include "calculate_integral.h"

#include <cmath>
#include <cstddef>
#include <vector>

#include "util/create_system_matrix.h"

namespace semester5_task4_1 {

double f(double x) {
    return std::sin(x);
}

std::pair<std::vector<std::vector<double>>, std::vector<double>> CreateSystem(
        std::vector<double> const& points) {
    std::vector<std::vector<double>> matrix = util::CreateSystemMatrix(points);

    std::vector<double> vector = std::vector<double>(points.size());
    for (std::size_t i = 0; i != vector.size(); ++i) {
        vector[i] = 1.0 / ((i + 2) * (i + 2));
    }
    return {matrix, vector};
}

}  // namespace semester5_task4_1
