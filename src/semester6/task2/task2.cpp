#include "task2.h"

#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

#include "model/matrix.h"
#include "solve_system.h"
#include "util/table.h"

namespace semester6_task2 {
void PrintSystem(model::Matrix const& matrix, std::vector<double> vector) {
    if (matrix.Size() <= 10) {
        std::cout << "Matrix:" << '\n';
        std::cout << matrix.ToString() << '\n';
    } else {
        std::cout << "Matrix of size " << matrix.Size() << ":" << '\n';
    }
    std::cout << "Condition number: " << matrix.NormConditionNumber() << '\n';

    auto const& [l, u] = matrix.LUDecomposition();
    if (matrix.Size() <= 10) {
        std::cout << "LU decomposition:" << '\n';
        std::cout << l->ToString() << u->ToString() << '\n';
    }
    std::cout << "L condition number: " << l->NormConditionNumber() << '\n';
    std::cout << "U condition number: " << u->NormConditionNumber() << '\n';

    auto const& [q, r] = matrix.QRDecomposition();
    if (matrix.Size() <= 10) {
        std::cout << "QR decomposition:" << '\n';
        std::cout << q->ToString() << r->ToString() << '\n';
    }
    std::cout << "Q condition number: " << q->NormConditionNumber() << '\n';
    std::cout << "R condition number: " << r->NormConditionNumber() << '\n';

    std::vector<std::vector<double>> table = SolveSystem(matrix, vector);
    if (matrix.Size() > 10) {
        std::cout << "First 10 rows of the table:" << '\n';
    }
    util::PrintTable(table, {"Right-hand side", "Library solution", "Gauss elimination",
                             "LU decomposition", "QR decomposition"});
    std::cout << '\n';
}

}  // namespace semester6_task2

namespace tasks {
void Semester6Task2() {
    std::cout << "Direct methods for solving systems of linear equations" << '\n';

    model::Matrix matrix1 = {{1, 7, 2}, {3, 5, 4}, {2, 8, 6}};
    std::vector<double> vector1 = {3, 7, 2};
    semester6_task2::PrintSystem(matrix1, std::move(vector1));

    model::Matrix matrix2 = {{3.278164, 1.045683, -1.378574},
                             {1.046583, 2.975937, 0.934251},
                             {-1.378574, 0.934251, 4.836173}};
    std::vector<double> vector2 = {-0.527466, 2.526877, 5.165441};
    matrix2.CacheLUDecomposition();
    matrix2.CacheQRDecomposition();
    semester6_task2::PrintSystem(matrix2, std::move(vector2));

    model::Matrix matrix3 = model::Matrix::CreateDiagonal(3, 2);
    matrix3.CacheLUDecomposition();
    matrix3.CacheQRDecomposition();
    semester6_task2::PrintSystem(matrix3);

    model::Matrix matrix4 = model::Matrix::CreateGilbert(3);
    matrix4.CacheLUDecomposition();
    matrix4.CacheQRDecomposition();
    semester6_task2::PrintSystem(matrix4);

    model::Matrix matrix5 = model::Matrix::CreateGilbert(6);
    matrix5.CacheLUDecomposition();
    matrix5.CacheQRDecomposition();
    semester6_task2::PrintSystem(matrix5);

    model::Matrix matrix6 = model::Matrix::CreateTridiagonal(4);
    matrix6.CacheLUDecomposition();
    matrix6.CacheQRDecomposition();
    semester6_task2::PrintSystem(matrix6);

    for (std::size_t i = 0; i != 5; ++i) {
        model::Matrix matrix7 = model::Matrix::CreateRandom(5);
        matrix7.CacheLUDecomposition();
        matrix7.CacheQRDecomposition();
        semester6_task2::PrintSystem(matrix7);
    }
}

}  // namespace tasks
