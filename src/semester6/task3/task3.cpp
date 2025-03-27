#include "task3.h"

#include <iostream>
#include <vector>

#include "model/matrix.h"
#include "solve_system.h"
#include "util/table.h"

namespace semester6_task3 {
void PrintSystem(model::Matrix const& matrix, double eps, std::vector<double> vector) {
    std::cout << "Matrix:" << '\n';
    std::cout << matrix.ToString() << '\n';

    std::vector<std::vector<double>> table = semester6_task3::SolveSystem(matrix, eps, vector);

    util::PrintTable(table, {"Right-hand side", "LU decomposition", "Jacobi iteration",
                             "Jacobi error", "Seidel iteration", "Seidel error"});
    std::cout << '\n';
}

}  // namespace semester6_task3

namespace tasks {
void Semester6Task3() {
    model::Matrix matrix1 = {{27, 7, 2}, {3, -30, 4}, {2, 8, 29}};
    std::vector<double> vector1 = {3, 7, 2};

    semester6_task3::PrintSystem(matrix1, 1e-5, vector1);
}

}  // namespace tasks
