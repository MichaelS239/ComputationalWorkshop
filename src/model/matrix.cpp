#include "matrix.h"

#include <cmath>
#include <cstddef>
#include <random>
#include <sstream>
#include <vector>

namespace model {
void Matrix::CalculateDeterminant(std::vector<bool> available_indices, unsigned depth,
                                  double& result) const {
    result = 0;
    int sign;
    std::size_t available_index = 0;
    for (std::size_t index = 0; index < matrix_.size(); index++) {
        if (available_indices[index]) {
            sign = (available_index % 2 == 0) ? 1 : -1;
            available_index++;
            available_indices[index] = false;
            double rec;
            CalculateDeterminant(available_indices, depth + 1, rec);
            if (depth != matrix_.size() - 1) {
                result += sign * matrix_[depth][index] * rec;
            } else
                result = sign * matrix_[depth][index];
            available_indices[index] = true;
        }
    }
}

double Matrix::Determinant() const {
    double det;
    CalculateDeterminant(std::vector<bool>(matrix_.size(), true), 0U, det);
    return det;
}

double Matrix::Norm() const {
    double norm = 0;
    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        for (std::size_t j = 0; j != matrix_.size(); ++j) {
            norm += matrix_[i][j] * matrix_[i][j];
        }
    }
    norm = std::sqrt(norm);
    return norm;
}

Matrix Matrix::Inverse() const {
    return util::CalculateInverseMatrix(matrix_);
}

std::vector<double> Matrix::SolveSystem(std::vector<double> const& vector) const {
    return util::SolveSystem(matrix_, vector);
}

double Matrix::NormConditionNumber() const {
    double norm = Norm();

    Matrix inverse_matrix = Inverse();
    double inverse_norm = inverse_matrix.Norm();

    return norm * inverse_norm;
}

double Matrix::VolumeConditionNumber() const {
    double det = Determinant();

    double product = 1;
    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        double cur_sum = 0;
        for (std::size_t j = 0; j != matrix_.size(); ++j) {
            cur_sum += matrix_[i][j] * matrix_[i][j];
        }
        product *= std::sqrt(cur_sum);
    }

    return product / std::abs(det);
}

double Matrix::AngleConditionNumber() const {
    double result = 0;

    Matrix inverse_matrix = util::CalculateInverseMatrix(matrix_);
    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        double scalar_product = 0;
        for (std::size_t j = 0; j != matrix_.size(); ++j) {
            scalar_product += std::abs(matrix_[i][j]) * std::abs(inverse_matrix[j][i]);
        }
        result = std::max(result, scalar_product);
    }

    return result;
}

std::string Matrix::ToString() const {
    std::stringstream test_ss;
    test_ss.precision(15);
    std::vector<std::vector<std::size_t>> sizes;
    std::vector<std::size_t> max_sizes(matrix_.size(), 0);
    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        sizes.emplace_back();
        for (std::size_t j = 0; j != matrix_[i].size(); ++j) {
            test_ss << matrix_[i][j];
            std::size_t size = test_ss.str().size();
            test_ss.str("");
            test_ss.clear();
            max_sizes[j] = std::max(max_sizes[j], size);
            sizes[sizes.size() - 1].push_back(size);
        }
    }
    std::stringstream ss;
    ss.precision(15);
    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        ss << '+';
        for (std::size_t j = 0; j != max_sizes[i] + 2; ++j) ss << '-';
    }
    ss << "+\n";
    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        for (std::size_t j = 0; j != matrix_[i].size(); ++j) {
            ss << "| " << matrix_[i][j];
            for (std::size_t k = sizes[i][j]; k != max_sizes[j] + 1; ++k) {
                ss << ' ';
            }
        }
        ss << "|\n";
    }

    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        ss << '+';
        for (std::size_t j = 0; j != max_sizes[i] + 2; ++j) ss << '-';
    }
    ss << "+\n";
    return ss.str();
}

Matrix Matrix::CreateDiagonal(std::size_t n, double k) {
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));
    for (std::size_t i = 0; i != n; ++i) {
        matrix[i][i] = k;
    }
    return matrix;
}

Matrix Matrix::CreateGilbert(std::size_t n) {
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));
    for (std::size_t i = 0; i != n; ++i) {
        for (std::size_t j = 0; j != n; ++j) {
            matrix[i][j] = 1 / (double)(i + j + 1);
        }
    }
    return matrix;
}

Matrix Matrix::CreateTridiagonal(std::size_t n) {
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));
    for (std::size_t i = 0; i != n; ++i) {
        for (std::size_t j = 0; j != n; ++j) {
            if (i == j) {
                matrix[i][j] = 2;
            } else if (i + 1 == j || i - 1 == j) {
                matrix[i][j] = -1;
            }
        }
    }
    return matrix;
}

Matrix Matrix::CreateRandom(std::size_t n, double lower_bound, double upper_bound) {
    std::vector<std::vector<double>> matrix(n, std::vector<double>(n));

    std::random_device r;
    std::default_random_engine re(r());
    std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    for (std::size_t i = 0; i != n; ++i) {
        for (std::size_t j = 0; j != n; ++j) {
            matrix[i][j] = unif(re);
        }
    }
    return matrix;
}

}  // namespace model
