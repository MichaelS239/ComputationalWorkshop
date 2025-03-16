#include "condition_numbers.h"

#include <cmath>
#include <cstddef>
#include <random>
#include <utility>
#include <vector>

#include "util/matrix_util.h"

namespace semester6_task1 {

std::pair<std::vector<double>, std::vector<std::vector<double>>> CalculateConds(
        model::Matrix const& matrix) {
    double norm_cond = matrix.NormConditionNumber();
    double volume_cond = matrix.VolumeConditionNumber();
    double angle_cond = matrix.AngleConditionNumber();

    std::vector<std::vector<double>> stats_table(matrix.Size(), std::vector<double>(6));
    std::vector<double> vec(matrix.Size());
    std::random_device r;
    std::default_random_engine re(r());
    std::uniform_real_distribution<double> unif(-10, 10);
    for (std::size_t i = 0; i != vec.size(); ++i) {
        vec[i] = unif(re);
        stats_table[i][0] = vec[i];
    }
    std::vector<double> first_solution = matrix.SolveSystem(vec);
    for (std::size_t i = 0; i != vec.size(); ++i) {
        stats_table[i][1] = first_solution[i];
    }
    std::random_device r1;
    std::default_random_engine re1(r());
    std::uniform_real_distribution<double> unif1(-0.1, 0.1);
    std::vector<double> new_vec(matrix.Size());
    for (std::size_t i = 0; i != vec.size(); ++i) {
        new_vec[i] = vec[i] + unif1(re1);
        stats_table[i][2] = new_vec[i];
    }
    std::vector<double> new_solution = matrix.SolveSystem(new_vec);
    for (std::size_t i = 0; i != vec.size(); ++i) {
        stats_table[i][3] = new_solution[i];
    }
    for (std::size_t i = 0; i != vec.size(); ++i) {
        double dif = std::abs(new_vec[i] - vec[i]);
        stats_table[i][4] = dif;
    }
    for (std::size_t i = 0; i != vec.size(); ++i) {
        double dif = std::abs(new_solution[i] - first_solution[i]);
        stats_table[i][5] = dif;
    }

    return {{norm_cond, volume_cond, angle_cond}, std::move(stats_table)};
}

}  // namespace semester6_task1
