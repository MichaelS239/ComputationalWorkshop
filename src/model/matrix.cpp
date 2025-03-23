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

double Matrix::Determinant(bool use_lu) const {
    if (!use_lu) {
        double det;
        CalculateDeterminant(std::vector<bool>(matrix_.size(), true), 0U, det);
        return det;
    }

    auto const& [l, u] = LUDecomposition();
    double det = 1;
    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        det *= (*l.get())[i][i];
    }
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

Matrix Matrix::Inverse(bool use_lu) const {
    if (!use_lu) return util::CalculateInverseMatrix(matrix_);

    Matrix inverse{matrix_.size()};
    std::vector<double> vector(matrix_.size());

    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        auto const& [l, u] = LUDecomposition();
        if (i != 0) vector[i - 1] = 0;
        vector[i] = 1;
        std::vector<double> solution = SolveSystem(vector);
        for (std::size_t j = 0; j != matrix_.size(); ++j) {
            inverse[j][i] = solution[j];
        }
    }

    return inverse;
}

std::vector<double> Matrix::SolveSystem(std::vector<double> const& vector,
                                        SolveMethod const solve_method) const {
    std::vector<double> solution;
    switch (solve_method) {
        case model::SolveMethod::Library: {
            solution = util::SolveSystem(matrix_, vector);
            break;
        }
        case model::SolveMethod::GaussElimination: {
            auto const& [new_matrix, new_vector] = GaussElimination(vector);
            solution = new_matrix.SolveUpperTriangularSystem(new_vector);
            break;
        }
        case model::SolveMethod::LUDecomposition: {
            auto const& [l, u] = LUDecomposition();
            std::vector<double> first_solution = l->SolveLowerTriangularSystem(vector);
            solution = u->SolveUpperTriangularSystem(first_solution);
            break;
        }
        case model::SolveMethod::QRDecomposition:
            break;
        default:
            break;
    }

    return solution;
}

std::vector<double> Matrix::SolveUpperTriangularSystem(std::vector<double> const& vector) const {
    std::vector<double> solution = vector;
    for (std::size_t k = matrix_.size() - 1; k != 0; --k) {
        solution[k] /= matrix_[k][k];
        for (std::size_t i = 0; i != k; ++i) {
            solution[i] -= solution[k] * matrix_[i][k];
        }
    }
    return solution;
}

std::vector<double> Matrix::SolveLowerTriangularSystem(std::vector<double> const& vector) const {
    std::vector<double> solution = vector;
    for (std::size_t k = 0; k != matrix_.size(); ++k) {
        solution[k] /= matrix_[k][k];
        for (std::size_t i = k + 1; i != matrix_.size(); ++i) {
            solution[i] -= solution[k] * matrix_[i][k];
        }
    }
    return solution;
}

std::pair<Matrix, std::vector<double>> Matrix::GaussElimination(
        std::vector<double> const& vector) const {
    Matrix matrix_copy{matrix_};
    std::vector<double> vector_copy = vector;

    for (std::size_t k = 0; k != matrix_.size(); ++k) {
        double tmp = matrix_copy[k][k];
        for (std::size_t i = k; i != matrix_.size(); ++i) {
            matrix_copy[k][i] /= tmp;
        }
        vector_copy[k] /= tmp;
        for (std::size_t i = k + 1; i != matrix_.size(); ++i) {
            tmp = matrix_copy[i][k];
            for (std::size_t j = k; j != matrix_.size(); ++j) {
                matrix_copy[i][j] -= matrix_copy[k][j] * tmp;
            }
            vector_copy[i] -= vector_copy[k] * tmp;
        }
    }

    return {std::move(matrix_copy), std::move(vector_copy)};
}

std::pair<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>> Matrix::LUDecomposition() const {
    if (l_cache_ != nullptr && u_cache_ != nullptr) return {l_cache_, u_cache_};

    Matrix l{matrix_.size()};
    Matrix u{matrix_.size()};

    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        l[i][0] = matrix_[i][0];
        u[0][i] = matrix_[0][i] / l[0][0];
    }

    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        for (std::size_t j = i; j != matrix_.size(); ++j) {
            l[j][i] = matrix_[j][i];
            for (std::size_t k = 0; k != i; ++k) {
                l[j][i] -= l[j][k] * u[k][i];
            }
        }
        for (std::size_t j = i; j != matrix_.size(); ++j) {
            u[i][j] = matrix_[i][j];
            for (std::size_t k = 0; k != i; ++k) {
                u[i][j] -= l[i][k] * u[k][j];
            }
            u[i][j] /= l[i][i];
        }
    }

    return {std::make_shared<Matrix>(l), std::make_shared<Matrix>(u)};
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
