#include "task5.h"

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <vector>

#include "model/matrix.h"

namespace semester6_task5 {
void PrintEigen(model::Matrix const& matrix, std::vector<double> const& eps) {
    if (matrix.Size() <= 10) {
        std::cout << "Matrix:" << '\n';
        std::cout << matrix.ToString() << '\n';
    } else {
        std::cout << "Matrix of size " << matrix.Size() << ":" << '\n';
    }

    for (std::size_t k = 0; k != eps.size(); ++k) {
        std::cout << "Accuracy: " << eps[k] << '\n';

        auto max_value_eigenvalues = matrix.GetEigenValues(eps[k]);
        std::sort(max_value_eigenvalues.eigenvalues.begin(),
                  max_value_eigenvalues.eigenvalues.end());
        if (matrix.Size() <= 10) {
            std::cout << "Eigenvalues of max value elimination: ";
            for (std::size_t i = 0; i != max_value_eigenvalues.eigenvalues.size(); ++i) {
                std::cout << std::setprecision(15) << max_value_eigenvalues.eigenvalues[i] << " ";
            }
            std::cout << '\n';
        } else {
            std::cout << "First five eigenvalues of max value elimination: ";
            for (std::size_t i = 0; i != 5; ++i) {
                std::cout << std::setprecision(15) << max_value_eigenvalues.eigenvalues[i] << " ";
            }
            std::cout << '\n';
            std::cout << "Last five eigenvalues of max value elimination: ";
            for (std::size_t i = max_value_eigenvalues.eigenvalues.size() - 5;
                 i != max_value_eigenvalues.eigenvalues.size(); ++i) {
                std::cout << std::setprecision(15) << max_value_eigenvalues.eigenvalues[i] << " ";
            }
            std::cout << '\n';
        }
        std::cout << "Number of iterations of max value elimination: "
                  << max_value_eigenvalues.iteration_num << '\n';

        auto cyclic_eigenvalues =
                matrix.GetEigenValues(eps[k], model::EliminationChoiceMethod::Cyclic);
        std::sort(cyclic_eigenvalues.eigenvalues.begin(), cyclic_eigenvalues.eigenvalues.end());
        if (matrix.Size() <= 10) {
            std::cout << "Eigenvalues of cyclic elimination: ";
            for (std::size_t i = 0; i != cyclic_eigenvalues.eigenvalues.size(); ++i) {
                std::cout << std::setprecision(15) << cyclic_eigenvalues.eigenvalues[i] << " ";
            }
            std::cout << '\n';
        } else {
            std::cout << "First five eigenvalues of cyclic elimination: ";
            for (std::size_t i = 0; i != 5; ++i) {
                std::cout << std::setprecision(15) << cyclic_eigenvalues.eigenvalues[i] << " ";
            }
            std::cout << '\n';
            std::cout << "Last five eigenvalues of cyclic elimination: ";
            for (std::size_t i = cyclic_eigenvalues.eigenvalues.size() - 5;
                 i != cyclic_eigenvalues.eigenvalues.size(); ++i) {
                std::cout << std::setprecision(15) << cyclic_eigenvalues.eigenvalues[i] << " ";
            }
            std::cout << '\n';
        }
        std::cout << "Number of iterations of cyclic elimination: "
                  << cyclic_eigenvalues.iteration_num << '\n';

        auto const& max_eigenvalue = matrix.MaxAbsoluteEigenvalue(eps[k]);
        auto const& min_eigenvalue = matrix.MinAbsoluteEigenvalue(eps[k]);

        std::cout << "Maximal absolute eigenvalue: " << std::setprecision(15)
                  << max_eigenvalue.eigenvalue << '\n';
        std::cout << "Minimal absolute eigenvalue: " << std::setprecision(15)
                  << min_eigenvalue.eigenvalue << '\n';
        std::cout << '\n';
    }
}

}  // namespace semester6_task5

namespace tasks {
void Semester6Task5() {
    model::Matrix matrix1 = {{-0.81417, -0.01937, 0.41372},
                             {-0.01937, 0.54414, 0.00590},
                             {0.41372, 0.00590, -0.81445}};
    semester6_task5::PrintEigen(matrix1, {1e-5, 1e-10});

    model::Matrix matrix2 = {{0.22, 0.02, 0.12, 0.14},
                             {0.02, 0.14, 0.04, -0.06},
                             {0.12, 0.04, 0.28, 0.08},
                             {0.14, -0.06, 0.08, 0.26}};
    semester6_task5::PrintEigen(matrix2, {1e-2, 1e-5, 1e-10});

    model::Matrix matrix3 = model::Matrix::CreateGilbert(5);
    semester6_task5::PrintEigen(matrix3, {1e-5, 1e-10});

    model::Matrix matrix4 = model::Matrix::CreateRandomSymmetric(100);
    semester6_task5::PrintEigen(matrix4, {1e-5, 1e-7});

    model::Matrix matrix5 = model::Matrix::CreateRandomSymmetric(300);
    semester6_task5::PrintEigen(matrix5, {1e-2, 1e-3});
}

}  // namespace tasks
