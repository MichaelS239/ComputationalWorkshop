#include "solve_system.h"

#include <random>
#include <vector>

namespace semester6_task2 {

std::vector<std::vector<double>> SolveSystem(model::Matrix const& matrix,
                                             std::vector<double> vector) {
    std::vector<std::vector<double>> table(matrix.Size(), std::vector<double>(5));

    if (vector.empty()) {
        vector.resize(matrix.Size());
        std::random_device rand;
        std::default_random_engine re(rand());
        std::uniform_real_distribution<double> unif(-10, 10);
        for (std::size_t i = 0; i != vector.size(); ++i) {
            vector[i] = unif(re);
        }
    }

    std::vector<double> library_solution = matrix.SolveSystem(vector);
    std::vector<double> gauss_solution =
            matrix.SolveSystem(vector, model::SolveMethod::GaussElimination);
    std::vector<double> lu_solution =
            matrix.SolveSystem(vector, model::SolveMethod::LUDecomposition);
    std::vector<double> qr_solution =
            matrix.SolveSystem(vector, model::SolveMethod::QRDecomposition);
    for (std::size_t i = 0; i != vector.size(); ++i) {
        table[i] = {vector[i], library_solution[i], gauss_solution[i], lu_solution[i],
                    qr_solution[i]};
    }

    return table;
}

}  // namespace semester6_task2
