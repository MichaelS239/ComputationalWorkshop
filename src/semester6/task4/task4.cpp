#include "task4.h"

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <vector>

#include "calculate_eigenvector.h"
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

        EigenInfo eigen_info = CalculateEigenvector(matrix, eps[k]);

        std::cout << "Power method:" << '\n';
        std::cout << '\n';
        std::cout << "Number of iterations of power method: " << eigen_info.power_iter_num << '\n';
        std::cout << "Maximal absolute eigenvalue of power method: " << std::setprecision(15)
                  << eigen_info.power_eigenvalue << '\n';
        std::cout << "Minimal absolute eigenvalue of power method: " << std::setprecision(15)
                  << eigen_info.power_min_eigenvalue << '\n';
        std::cout << "Check of eigenvector of power method:" << '\n';
        util::PrintTable(eigen_info.power_table, {"Eigenvector (x)", "A * x", "(A * x) / x"});
        std::cout << '\n';
        std::cout << "Scalar product method:" << '\n';
        std::cout << "Number of iterations of scalar product method: " << eigen_info.scalar_iter_num
                  << '\n';
        std::cout << "Maximal absolute eigenvalue of scalar product method: "
                  << std::setprecision(15) << eigen_info.scalar_eigenvalue << '\n';
        std::cout << "Minimal absolute eigenvalue of scalar product method: "
                  << std::setprecision(15) << eigen_info.scalar_min_eigenvalue << '\n';
        std::cout << "Check of eigenvector of scalar product method:" << '\n';
        util::PrintTable(eigen_info.scalar_table, {"Eigenvector (x)", "A * x", "(A * x) / x"});
        std::cout << '\n';
    }
}

}  // namespace semester6_task4

namespace tasks {
void Semester6Task4() {
    std::cout << "Numerical methods for finding the maximal absolute eigenvalue of a matrix"
              << '\n';
    std::cout << '\n';
    model::Matrix matrix1 = {{-0.81417, -0.01937, 0.41372},
                             {-0.01937, 0.54414, 0.00590},
                             {0.41372, 0.00590, -0.81445}};
    semester6_task4::PrintEigen(matrix1, {1e-5, 1e-10});

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

    model::Matrix matrix5 = model::Matrix::CreateGilbert(5);
    semester6_task4::PrintEigen(matrix5, {1e-5, 1e-10});
}

}  // namespace tasks
