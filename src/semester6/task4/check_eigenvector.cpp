#include "check_eigenvector.h"

namespace semester6_task4 {
std::vector<std::vector<double>> CheckEigenvector(model::Matrix const& matrix,
                                                  std::vector<double> const& eigenvector) {
    std::vector<std::vector<double>> table(matrix.Size(), std::vector<double>(4));

    std::vector<double> multiplied_eigenvector = matrix.MultiplyByVector(eigenvector);
    for (std::size_t i = 0; i != matrix.Size(); ++i) {
        table[i] = {eigenvector[i], multiplied_eigenvector[i],
                    multiplied_eigenvector[i] / eigenvector[i]};
    }
    return table;
}

}  // namespace semester6_task4
