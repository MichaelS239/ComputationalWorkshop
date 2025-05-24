#include "task1.h"

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "condition_numbers.h"
#include "model/matrix.h"
#include "util/table.h"

namespace semester6_task1 {
void PrintConds(model::Matrix const &matrix) {
    std::cout << "Matrix:" << '\n';
    std::cout << matrix.ToString() << '\n';
    auto const &[conds, stats_table] = semester6_task1::CalculateConds(matrix);

    std::cout << "Norm condition number: " << conds[0] << '\n';
    std::cout << "Volume condition number: " << conds[1] << '\n';
    std::cout << "Angle condition number: " << conds[2] << '\n';
    std::cout << '\n';

    util::PrintTable(stats_table,
                     {"Right-hand side", "Solution", "New right-hand side", "New solution",
                      "Right-hand side difference", "Solution difference"});
    std::cout << '\n';
}

}  // namespace semester6_task1

namespace tasks {
void Semester6Task1() {
    std::cout << "Condition numbers of a matrix" << '\n';
    std::cout << '\n';
    model::Matrix matrix1 = model::Matrix::CreateDiagonal(3, 2);
    semester6_task1::PrintConds(matrix1);
    model::Matrix matrix2 = model::Matrix::CreateDiagonal(4, 2);
    semester6_task1::PrintConds(matrix2);
    model::Matrix matrix3 = model::Matrix::CreateDiagonal(10, 2);
    semester6_task1::PrintConds(matrix3);
    model::Matrix matrix4 = model::Matrix::CreateGilbert(3);
    semester6_task1::PrintConds(matrix4);
    model::Matrix matrix5 = model::Matrix::CreateGilbert(4);
    semester6_task1::PrintConds(matrix5);
    model::Matrix matrix6 = model::Matrix::CreateGilbert(8);
    semester6_task1::PrintConds(matrix6);
    model::Matrix matrix7 = model::Matrix::CreateTridiagonal(4);
    semester6_task1::PrintConds(matrix7);
    model::Matrix matrix8 = model::Matrix::CreateTridiagonal(8);
    semester6_task1::PrintConds(matrix8);
    model::Matrix matrix9 = model::Matrix::CreateTridiagonal(10);
    semester6_task1::PrintConds(matrix9);
    model::Matrix matrix10 = model::Matrix::CreateRandom(5);
    semester6_task1::PrintConds(matrix10);
}

}  // namespace tasks
