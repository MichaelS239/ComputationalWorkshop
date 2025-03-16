#include "matrix_util.h"

#include <Eigen/Dense>
#include <cstddef>
#include <vector>

namespace util {

using Eigen::MatrixXd;
using Eigen::VectorXd;

std::vector<double> SolveSystem(std::vector<std::vector<double>> const& matrix,
                                std::vector<double> const& vector) {
    MatrixXd m(matrix.size(), matrix[0].size());
    for (std::size_t i = 0; i != matrix.size(); ++i) {
        for (std::size_t j = 0; j != matrix.size(); ++j) {
            m(i, j) = matrix[i][j];
        }
    }
    VectorXd v(vector.size());
    for (std::size_t i = 0; i != vector.size(); ++i) {
        v(i) = vector[i];
    }
    VectorXd s = m.colPivHouseholderQr().solve(v);
    std::vector<double> solution = std::vector<double>(matrix.size());
    for (std::size_t i = 0; i != matrix.size(); ++i) {
        solution[i] = s(i);
    }
    return solution;
}

std::vector<std::vector<double>> CalculateInverseMatrix(
        std::vector<std::vector<double>> const& matrix) {
    MatrixXd m(matrix.size(), matrix[0].size());
    for (std::size_t i = 0; i != matrix.size(); ++i) {
        for (std::size_t j = 0; j != matrix.size(); ++j) {
            m(i, j) = matrix[i][j];
        }
    }
    MatrixXd inverse = m.inverse();
    std::vector<std::vector<double>> inverse_matrix(matrix.size(),
                                                    std::vector<double>(matrix[0].size()));
    for (std::size_t i = 0; i != matrix.size(); ++i) {
        for (std::size_t j = 0; j != matrix.size(); ++j) {
            inverse_matrix[i][j] = inverse(i, j);
        }
    }
    return inverse_matrix;
}

}  // namespace util
