#include "create_system_matrix.h"

#include <cstddef>
#include <vector>

namespace util {
std::vector<std::vector<double>> CreateSystemMatrix(std::vector<double> const& points) {
    std::vector<std::vector<double>> matrix =
            std::vector<std::vector<double>>(points.size(), std::vector<double>(points.size()));
    for (std::size_t i = 0; i != matrix.size(); ++i) {
        for (std::size_t j = 0; j != matrix.size(); ++j) {
            if (i == 0) {
                matrix[i][j] = 1;
            } else {
                matrix[i][j] = matrix[i - 1][j] * points[j];
            }
        }
    }

    return matrix;
}

}  // namespace util
