#include "solve_system.h"

#include <random>
#include <vector>

namespace semester6_task3 {

std::vector<std::vector<double>> SolveSystem(model::Matrix const& matrix, double eps,
                                             std::vector<double> vector) {
    std::vector<std::vector<double>> table(matrix.Size(), std::vector<double>(4));

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
    std::vector<double> jacobi_solution =
            matrix.SolveSystem(vector, eps, model::IterationMethod::Jacobi);
    for (std::size_t i = 0; i != vector.size(); ++i) {
        table[i] = {vector[i], lu_solution[i], jacobi_solution[i],
                    std::abs(lu_solution[i] - jacobi_solution[i])};
    }

    return table;
}

}  // namespace semester6_task3
