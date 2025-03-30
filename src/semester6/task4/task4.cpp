#include "task4.h"

#include <cstddef>
#include <iostream>
#include <vector>

#include "check_eigenvector.h"
#include "model/matrix.h"
#include "util/table.h"

namespace semester6_task4 {
void PrintEigen(model::Matrix const& matrix, std::vector<double> const& eps) {
    if (matrix.Size() <= 10) {
        std::cout << "Matrix:" << '\n';
        std::cout << matrix.ToString() << '\n';
    } else {
        std::cout << "Matrix of size " << matrix.Size() << ":" << '\n';
    }

    for (std::size_t k = 0; k != eps.size(); ++k) {
        std::cout << "Accuracy: " << eps[k] << '\n';

        auto const& [eigen_pair, iter_num] = matrix.MaxEigenvalue(eps[k]);
        std::cout << "Number of iterations: " << iter_num << '\n';
        std::cout << "Eigenvalue: " << eigen_pair.first << '\n';
        std::vector<std::vector<double>> table = CheckEigenvector(matrix, eigen_pair.second);
        util::PrintTable(table, {"Eigenvector (x)", "A * x", "(A * x) / x"});
    }
}

}  // namespace semester6_task4

namespace tasks {
void Semester6Task4() {
    model::Matrix matrix1 = {{-0.81417, -0.01937, 0.41372},
                             {-0.01937, 0.54414, 0.00590},
                             {0.41372, 0.00590, -0.81445}};
    semester6_task4::PrintEigen(matrix1, {1e-5, 1e-7});

    model::Matrix matrix2 = {{0.22, 0.02, 0.12, 0.14},
                             {0.02, 0.14, 0.04, -0.06},
                             {0.12, 0.04, 0.28, 0.08},
                             {0.14, -0.06, 0.08, 0.26}};
    semester6_task4::PrintEigen(matrix2, {1e-2, 1e-5, 1e-10});

    model::Matrix matrix3 = {{0.3, 1.4, -0.2, 0.5},
                             {-2.1, 1.2, 1.7, 1.2},
                             {3.2, 0.1, -1.4, 0.1},
                             {0.8, -0.2, 1.5, -0.1}};
    semester6_task4::PrintEigen(matrix3, {1e-2, 1e-5, 1e-10});

    model::Matrix matrix4 = {{1.3, 0.2, -0.1, 0.8},
                             {-0.4, 1.5, 0.2, 1.2},
                             {1.1, -1.8, 1.3, -0.2},
                             {0.8, -0.1, 0.2, 1.0}};
    semester6_task4::PrintEigen(matrix4, {1e-2, 1e-5, 1e-10});
}

}  // namespace tasks
