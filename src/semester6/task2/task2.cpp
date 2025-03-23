#include "task2.h"

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "model/matrix.h"
#include "util/table.h"

namespace tasks {
void Semester6Task2() {
    std::cout << "Direct methods for solving systems of linear equations" << '\n';
    model::Matrix matrix = {{1, 7, 2}, {3, 5, 4}, {2, 8, 6}};
    std::vector<double> vector = {3, 7, 2};
    std::vector<double> first_solution = matrix.SolveSystem(vector);
    std::vector<double> second_solution =
            matrix.SolveSystem(vector, model::SolveMethod::GaussElimination);
    std::vector<double> third_solution =
            matrix.SolveSystem(vector, model::SolveMethod::LUDecomposition);
    std::vector<double> forth_solution =
            matrix.SolveSystem(vector, model::SolveMethod::QRDecomposition);
    for (std::size_t i = 0; i != first_solution.size(); ++i) {
        std::cout << first_solution[i] << " ";
    }
    std::cout << '\n';
    for (std::size_t i = 0; i != second_solution.size(); ++i) {
        std::cout << second_solution[i] << " ";
    }
    std::cout << '\n';
    for (std::size_t i = 0; i != third_solution.size(); ++i) {
        std::cout << third_solution[i] << " ";
    }
    std::cout << '\n';
    for (std::size_t i = 0; i != forth_solution.size(); ++i) {
        std::cout << forth_solution[i] << " ";
    }
    std::cout << '\n';

    model::Matrix inv1 = matrix.Inverse();
    model::Matrix inv2 = matrix.Inverse(model::SolveMethod::GaussElimination);
    model::Matrix inv3 = matrix.Inverse(model::SolveMethod::LUDecomposition);
    model::Matrix inv4 = matrix.Inverse(model::SolveMethod::QRDecomposition);
    std::cout << inv1.ToString() << '\n';
    std::cout << inv2.ToString() << '\n';
    std::cout << inv3.ToString() << '\n';
    std::cout << inv4.ToString() << '\n';

    std::cout << matrix.Determinant() << '\n';
    std::cout << matrix.Determinant(model::CalculationMethod::LUDecomposition) << '\n';
    std::cout << matrix.Determinant(model::CalculationMethod::QRDecomposition) << '\n';
}

}  // namespace tasks
