#include "solve_system.h"

#include <cstddef>
#include <random>
#include <utility>
#include <vector>

namespace semester6_task3 {

std::pair<std::vector<std::vector<double>>, std::pair<std::size_t, std::size_t>> SolveSystem(
        model::Matrix const& matrix, double eps, std::vector<double> vector) {
    std::vector<std::vector<double>> table(std::min(10UL, matrix.Size()), std::vector<double>(6));

    if (vector.empty()) {
        vector.resize(matrix.Size());
        std::random_device rand;
        std::default_random_engine re(rand());
        std::uniform_real_distribution<double> unif(-10, 10);
        for (std::size_t i = 0; i != vector.size(); ++i) {
            vector[i] = unif(re);
        }
    }

    std::vector<double> lu_solution =
            matrix.SolveSystem(vector, model::SolveMethod::LUDecomposition);
    auto const& [jacobi_solution, jacobi_iter_num] =
            matrix.SolveSystem(vector, eps, model::IterationMethod::Jacobi);
    auto const& [seidel_solution, seidel_iter_num] =
            matrix.SolveSystem(vector, eps, model::IterationMethod::Seidel);
    for (std::size_t i = 0; i != std::min(10UL, vector.size()); ++i) {
        table[i] = {vector[i],          lu_solution[i],
                    jacobi_solution[i], std::abs(lu_solution[i] - jacobi_solution[i]),
                    seidel_solution[i], std::abs(lu_solution[i] - seidel_solution[i])};
    }

    return {std::move(table), {jacobi_iter_num, seidel_iter_num}};
}

}  // namespace semester6_task3
