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
    for (std::size_t i = 0; i != first_solution.size(); ++i) {
        std::cout << first_solution[i] << " ";
    }
    std::cout << '\n';
    for (std::size_t i = 0; i != second_solution.size(); ++i) {
        std::cout << second_solution[i] << " ";
    }
    std::cout << '\n';
}

}  // namespace tasks
