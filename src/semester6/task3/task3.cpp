#include "task3.h"

#include <cstddef>
#include <iostream>
#include <vector>

#include "model/matrix.h"
#include "solve_system.h"
#include "util/table.h"

namespace semester6_task3 {
void PrintSystem(model::Matrix const& matrix, std::vector<double> const& eps,
                 std::vector<double> vector) {
    if (matrix.Size() <= 10) {
        std::cout << "Matrix:" << '\n';
        std::cout << matrix.ToString() << '\n';
    } else {
        std::cout << "Matrix of size " << matrix.Size() << ":" << '\n';
    }
    for (std::size_t i = 0; i != eps.size(); ++i) {
        std::cout << "Accuracy: " << eps[i] << '\n';

        auto const& [table, iter_pair] = semester6_task3::SolveSystem(matrix, eps[i], vector);
        std::cout << "Number of iterations of Jacobi method: " << iter_pair.first << '\n';
        std::cout << "Number of iterations of Seidel method: " << iter_pair.second << '\n';
        if (matrix.Size() > 10) {
            std::cout << "First 10 rows of the table:" << '\n';
        }
        util::PrintTable(table, {"Right-hand side", "LU decomposition", "Jacobi iteration",
                                 "Jacobi error", "Seidel iteration", "Seidel error"});
        std::cout << '\n';
    }
}

}  // namespace semester6_task3

namespace tasks {
void Semester6Task3() {
    model::Matrix matrix1 = {{27, 7, 2}, {3, -30, 4}, {2, 8, 29}};
    std::vector<double> vector1 = {3, 7, 2};

    semester6_task3::PrintSystem(matrix1, {1e-5, 1e-15}, vector1);

    model::Matrix matrix2 = model::Matrix::CreateDiagonal(100, 200);
    for (std::size_t i = 0; i != matrix2.Size() - 1; ++i) {
        matrix2[i + 1][i] = -1;
        matrix2[matrix2.Size() - 2 - i][matrix2.Size() - 1 - i] = -1;
    }
    semester6_task3::PrintSystem(matrix2, {1e-5, 1e-15});

    model::Matrix matrix3 = model::Matrix::CreateRandomDiagonalDominant(5);
    semester6_task3::PrintSystem(matrix3, {1e-5});

    model::Matrix matrix4 = model::Matrix::CreateRandomSymmetricDiagonalDominant(5);
    semester6_task3::PrintSystem(matrix4, {1e-5});

    model::Matrix matrix5 = model::Matrix::CreateRandomSymmetricDiagonalDominant(1000);
    semester6_task3::PrintSystem(matrix5, {1e-5, 1e-10});

    model::Matrix matrix6 = model::Matrix::CreateRandomSymmetricDiagonalDominant(2000);
    semester6_task3::PrintSystem(matrix6, {1e-5, 1e-10});
}

}  // namespace tasks
