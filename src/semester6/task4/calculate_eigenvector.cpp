#include "calculate_eigenvector.h"

#include <cstddef>
#include <utility>
#include <vector>

namespace semester6_task4 {
EigenInfo CalculateEigenvector(model::Matrix const& matrix, double eps) {
    std::vector<std::vector<double>> power_table(matrix.Size(), std::vector<double>(3));
    std::vector<std::vector<double>> scalar_table(matrix.Size(), std::vector<double>(3));

    auto const& power_eigen_info = matrix.MaxAbsoluteEigenvalue(eps);
    auto const& power_min_eigen_info = matrix.MinAbsoluteEigenvalue(eps);
    std::vector<double> multiplied_power_eigenvector =
            matrix.MultiplyByVector(power_eigen_info.eigenvector);

    auto const& scalar_eigen_info =
            matrix.MaxAbsoluteEigenvalue(eps, model::EigenvalueMethod::ScalarProduct);
    auto const& scalar_min_eigen_info =
            matrix.MinAbsoluteEigenvalue(eps, model::EigenvalueMethod::ScalarProduct);
    std::vector<double> multiplied_scalar_eigenvector =
            matrix.MultiplyByVector(scalar_eigen_info.eigenvector);
    for (std::size_t i = 0; i != matrix.Size(); ++i) {
        power_table[i] = {power_eigen_info.eigenvector[i], multiplied_power_eigenvector[i],
                          multiplied_power_eigenvector[i] / power_eigen_info.eigenvector[i]};
        scalar_table[i] = {scalar_eigen_info.eigenvector[i], multiplied_scalar_eigenvector[i],
                           multiplied_scalar_eigenvector[i] / scalar_eigen_info.eigenvector[i]};
    }
    return {std::move(power_table),          std::move(scalar_table),
            power_eigen_info.eigenvalue,     scalar_eigen_info.eigenvalue,
            power_eigen_info.iteration_num,  scalar_eigen_info.iteration_num,
            power_min_eigen_info.eigenvalue, scalar_min_eigen_info.eigenvalue};
}

}  // namespace semester6_task4
