#include "matrix.h"

#include <cmath>
#include <cstddef>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "util/matrix_util.h"
#include "util/vector_util.h"

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

double Matrix::Determinant(CalculationMethod calc_method) const {
    double det;
    switch (calc_method) {
        case model::CalculationMethod::Library:
            CalculateDeterminant(std::vector<bool>(matrix_.size(), true), 0U, det);
            break;
        case model::CalculationMethod::LUDecomposition: {
            auto const& [l, u] = LUDecomposition();
            det = 1;
            for (std::size_t i = 0; i != matrix_.size(); ++i) {
                det *= (*l.get())[i][i];
            }
            break;
        }
        case model::CalculationMethod::QRDecomposition: {
            auto const& [q, r] = QRDecomposition();
            det = 1;
            for (std::size_t i = 0; i != matrix_.size(); ++i) {
                det *= (*r.get())[i][i];
            }
            break;
        }
        default:
            break;
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

Matrix Matrix::Inverse(SolveMethod calc_method) const {
    if (calc_method == model::SolveMethod::Library) {
        return util::CalculateInverseMatrix(matrix_);
    }

    Matrix inverse{matrix_.size()};
    std::vector<double> vector(matrix_.size());

    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        if (i != 0) vector[i - 1] = 0;
        vector[i] = 1;
        std::vector<double> solution = SolveSystem(vector, calc_method);
        for (std::size_t j = 0; j != matrix_.size(); ++j) {
            inverse[j][i] = solution[j];
        }
    }

    return inverse;
}

Matrix Matrix::Transpose() const {
    Matrix matrix{matrix_};
    matrix.SelfTranspose();
    return matrix;
}

void Matrix::SelfTranspose() {
    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        for (std::size_t j = 0; j != i; ++j) {
            std::swap(matrix_[i][j], matrix_[j][i]);
        }
    }
}

std::vector<double> Matrix::MultiplyByVector(std::vector<double> const& vector) const {
    if (vector.size() != matrix_.size())
        throw std::runtime_error("Error: sizes of matrix and vector are different");

    std::vector<double> result(vector.size());
    for (std::size_t i = 0; i != vector.size(); ++i) {
        for (std::size_t k = 0; k != vector.size(); ++k) {
            result[i] += matrix_[i][k] * vector[k];
        }
    }

    return result;
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
        case model::SolveMethod::QRDecomposition: {
            auto const& [q, r] = QRDecomposition();
            std::vector<double> new_vector = (q->Transpose()).MultiplyByVector(vector);
            solution = r->SolveUpperTriangularSystem(new_vector);
            break;
        }
        default:
            break;
    }

    return solution;
}

std::vector<double> Matrix::SolveUpperTriangularSystem(std::vector<double> const& vector) const {
    std::vector<double> solution = vector;
    for (std::size_t k = matrix_.size() - 1; k >= 0; --k) {
        solution[k] /= matrix_[k][k];
        for (std::size_t i = 0; i != k; ++i) {
            solution[i] -= solution[k] * matrix_[i][k];
        }
        if (k == 0) break;
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

std::vector<double> Matrix::SolveTridiagonalSystem(std::vector<double> const& vector) const {
    std::vector<double> first_coefs(vector.size() + 1);
    std::vector<double> second_coefs(vector.size() + 1);
    std::vector<double> solution(vector.size());

    for (std::size_t i = 0; i != vector.size(); ++i) {
        if (i == 0) {
            first_coefs[i + 1] = -matrix_[i][i + 1] / matrix_[i][i];
            second_coefs[i + 1] = vector[i] / matrix_[i][i];
        } else if (i == vector.size() - 1) {
            second_coefs[i + 1] = (matrix_[i][i - 1] * second_coefs[i] - vector[i]) /
                                  (-matrix_[i][i] - matrix_[i][i - 1] * first_coefs[i]);
        } else {
            first_coefs[i + 1] =
                    matrix_[i][i + 1] / (-matrix_[i][i] - matrix_[i][i - 1] * first_coefs[i]);
            second_coefs[i + 1] = (matrix_[i][i - 1] * second_coefs[i] - vector[i]) /
                                  (-matrix_[i][i] - matrix_[i][i - 1] * first_coefs[i]);
        }
    }

    solution[vector.size() - 1] = second_coefs[vector.size()];
    for (std::size_t i = 0; i != vector.size() - 1; ++i) {
        solution[vector.size() - 2 - i] =
                first_coefs[vector.size() - 1 - i] * solution[vector.size() - 1 - i] +
                second_coefs[vector.size() - 1 - i];
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
        for (std::size_t j = i; j != matrix_.size(); ++j) {
            u[i][j] = matrix_[i][j];
            for (std::size_t k = 0; k != i; ++k) {
                u[i][j] -= l[i][k] * u[k][j];
            }
        }
        for (std::size_t j = i; j != matrix_.size(); ++j) {
            l[j][i] = matrix_[j][i];
            for (std::size_t k = 0; k != i; ++k) {
                l[j][i] -= l[j][k] * u[k][i];
            }
            l[j][i] /= u[i][i];
        }
    }

    return {std::make_shared<Matrix>(l), std::make_shared<Matrix>(u)};
}

void Matrix::MultiplyByOrthogonalLeft(std::size_t i, std::size_t j, double cosine, double sine) {
    for (std::size_t k = 0; k != matrix_.size(); ++k) {
        double new_i = matrix_[i][k] * cosine - matrix_[j][k] * sine;
        double new_j = matrix_[i][k] * sine + matrix_[j][k] * cosine;
        matrix_[i][k] = new_i;
        matrix_[j][k] = new_j;
    }
}

void Matrix::MultiplyByOrthogonalRight(std::size_t i, std::size_t j, double cosine, double sine) {
    for (std::size_t k = 0; k != matrix_.size(); ++k) {
        double new_i = matrix_[k][i] * cosine + matrix_[k][j] * sine;
        double new_j = -matrix_[k][i] * sine + matrix_[k][j] * cosine;
        matrix_[k][i] = new_i;
        matrix_[k][j] = new_j;
    }
}

std::pair<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>> Matrix::QRDecomposition() const {
    if (q_cache_ != nullptr && r_cache_ != nullptr) return {q_cache_, r_cache_};

    Matrix q = CreateDiagonal(matrix_.size(), 1);
    Matrix r{matrix_};

    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        for (std::size_t j = i + 1; j != matrix_.size(); ++j) {
            if (r[i][i] != 0 || r[j][i] != 0) {
                double root = std::sqrt(r[i][i] * r[i][i] + r[j][i] * r[j][i]);
                double cosine = r[i][i] / root;
                double sine = -r[j][i] / root;
                r.MultiplyByOrthogonalLeft(i, j, cosine, sine);
                q.MultiplyByOrthogonalLeft(i, j, cosine, sine);
            }
        }
    }
    q.SelfTranspose();

    return {std::make_shared<Matrix>(q), std::make_shared<Matrix>(r)};
}

void Matrix::CacheLUDecomposition() {
    auto const& [l, u] = LUDecomposition();
    l_cache_ = l;
    u_cache_ = u;
}

void Matrix::CacheQRDecomposition() {
    auto const& [q, r] = QRDecomposition();
    q_cache_ = q;
    r_cache_ = r;
}

std::pair<std::vector<double>, std::size_t> Matrix::SolveSystem(
        std::vector<double> const& vector, double eps, IterationMethod const solve_method) const {
    std::vector<double> solution;
    std::size_t iter_num = 0;
    switch (solve_method) {
        case model::IterationMethod::Jacobi: {
            auto [sol, iter] = JacobiIteration(vector, eps);
            solution = std::move(sol);
            iter_num = iter;
            break;
        }
        case model::IterationMethod::Seidel: {
            auto [sol, iter] = SeidelIteration(vector, eps);
            solution = std::move(sol);
            iter_num = iter;
            break;
        }
        default:
            break;
    }

    return {std::move(solution), iter_num};
}

std::pair<std::vector<double>, std::size_t> Matrix::JacobiIteration(
        std::vector<double> const& vector, double eps) const {
    if (!IsDiagonallyDominant()) {
        throw std::runtime_error("Error: matrix must be diagonally dominant");
    }

    Matrix b{matrix_.size()};
    std::vector<double> c(matrix_.size());
    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        for (std::size_t j = 0; j != matrix_.size(); ++j) {
            if (j != i) {
                b[i][j] = -matrix_[i][j] / matrix_[i][i];
            }
        }
        c[i] = vector[i] / matrix_[i][i];
    }

    std::vector<double> x0(matrix_.size());
    std::vector<double> x1(matrix_.size());
    double difference = 0;
    std::size_t iteration_num = 0;

    do {
        ++iteration_num;
        x0 = std::move(x1);
        difference = 0;
        x1 = b.MultiplyByVector(x0);
        for (std::size_t i = 0; i != x1.size(); ++i) {
            x1[i] += c[i];
            difference += (x1[i] - x0[i]) * (x1[i] - x0[i]);
        }
    } while (std::sqrt(difference) >= eps);

    return {std::move(x1), iteration_num};
}

std::pair<std::vector<double>, std::size_t> Matrix::SeidelIteration(
        std::vector<double> const& vector, double eps) const {
    std::vector<double> x0(matrix_.size());
    std::vector<double> x1(matrix_.size());
    double difference = 0;
    std::size_t iteration_num = 0;

    do {
        ++iteration_num;
        x0 = std::move(x1);
        difference = 0;
        x1 = std::vector<double>(matrix_.size());
        for (std::size_t i = 0; i != x1.size(); ++i) {
            for (std::size_t j = 0; j != i; ++j) {
                x1[i] -= matrix_[i][j] / matrix_[i][i] * x1[j];
            }
            for (std::size_t j = i + 1; j != matrix_.size(); ++j) {
                x1[i] -= matrix_[i][j] / matrix_[i][i] * x0[j];
            }
            x1[i] += vector[i] / matrix_[i][i];
            difference += (x1[i] - x0[i]) * (x1[i] - x0[i]);
        }
    } while (std::sqrt(difference) >= eps);

    return {std::move(x1), iteration_num};
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

Matrix::EigenInfo Matrix::PowerMethod(double eps) const {
    std::vector<double> x0(matrix_.size(), 1);
    std::vector<double> x1(matrix_.size(), 1);

    double diff = 0;
    double eigenvalue = 0;
    double old_eigenvalue = 0;
    std::size_t iteration_num = 0;
    do {
        ++iteration_num;
        x0 = std::move(x1);
        old_eigenvalue = eigenvalue;
        x1 = MultiplyByVector(x0);
        eigenvalue = util::ScalarProduct(x1, x0) / util::ScalarProduct(x0, x0);
        double diff_numerator = 0;
        for (std::size_t i = 0; i != matrix_.size(); ++i) {
            diff_numerator += (x1[i] - eigenvalue * x0[i]) * (x1[i] - eigenvalue * x0[i]);
        }
        diff = std::sqrt(diff_numerator) / util::Norm(x0);
        if (iteration_num % 10 == 0) {
            util::Normalize(x1);
        }
    } while (diff >= eps);

    util::Normalize(x1);

    return {eigenvalue, std::move(x1), iteration_num};
}

Matrix::EigenInfo Matrix::ScalarProductMethod(double eps) const {
    model::Matrix transposed_matrix = Transpose();
    std::vector<double> x0(matrix_.size(), 1);
    std::vector<double> x1(matrix_.size(), 1);
    std::vector<double> y0(matrix_.size(), 1);
    std::vector<double> y1(matrix_.size(), 1);

    double diff = 0;
    double eigenvalue = 0;
    double old_eigenvalue = 0;
    std::size_t iteration_num = 0;
    do {
        ++iteration_num;
        x0 = std::move(x1);
        y0 = std::move(y1);
        old_eigenvalue = eigenvalue;
        x1 = MultiplyByVector(x0);
        y1 = transposed_matrix.MultiplyByVector(y0);
        eigenvalue = util::ScalarProduct(x1, y1) / util::ScalarProduct(x0, y1);
        double diff_numerator = 0;
        for (std::size_t i = 0; i != matrix_.size(); ++i) {
            diff_numerator += (x1[i] - eigenvalue * x0[i]) * (x1[i] - eigenvalue * x0[i]);
        }
        diff = std::sqrt(diff_numerator) / util::Norm(x0);
        if (iteration_num % 10 == 0) {
            util::Normalize(x1);
            util::Normalize(y1);
        }
    } while (diff >= eps);

    util::Normalize(x1);

    return {eigenvalue, std::move(x1), iteration_num};
}

Matrix::EigenInfo Matrix::MaxAbsoluteEigenvalue(double eps,
                                                EigenvalueMethod const eigenvalue_method) const {
    EigenInfo eigen_info;
    switch (eigenvalue_method) {
        case model::EigenvalueMethod::Power: {
            auto eigen = PowerMethod(eps);
            eigen_info = std::move(eigen);
            break;
        }
        case model::EigenvalueMethod::ScalarProduct: {
            auto eigen = ScalarProductMethod(eps);
            eigen_info = std::move(eigen);
            break;
        }
        default:
            break;
    }

    return eigen_info;
}

Matrix::EigenInfo Matrix::MinAbsoluteEigenvalue(double eps,
                                                EigenvalueMethod const eigenvalue_method) const {
    Matrix inv = Inverse();
    Matrix::EigenInfo inv_info = inv.MaxAbsoluteEigenvalue(eps, eigenvalue_method);
    inv_info.eigenvalue = 1 / inv_info.eigenvalue;
    return inv_info;
}

double Matrix::CalculateNonDiagonalSum() const {
    double non_diagonal_sum = 0;
    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        for (std::size_t j = 0; j != matrix_.size(); ++j) {
            if (i != j) {
                non_diagonal_sum += matrix_[i][j] * matrix_[i][j];
            }
        }
    }

    return non_diagonal_sum;
}

Matrix::EigenValuesInfo Matrix::GetEigenValues(
        double eps, model::EliminationChoiceMethod const elimination_method) const {
    if (!IsSymmetric()) {
        throw std::runtime_error("Error: matrix must be symmetric");
    }

    Matrix matrix_copy{matrix_};

    double non_diagonal_sum = CalculateNonDiagonalSum();
    std::size_t iteration_num = 0;
    std::size_t max_i = 0;
    std::size_t max_j = 0;
    do {
        ++iteration_num;
        double max_abs = 0;

        switch (elimination_method) {
            case model::EliminationChoiceMethod::MaxValue: {
                max_i = 0;
                max_j = 0;
                for (std::size_t i = 0; i != matrix_copy.Size(); ++i) {
                    for (std::size_t j = 0; j != matrix_copy.Size(); ++j) {
                        if (i != j && std::abs(matrix_copy[i][j]) > max_abs) {
                            max_abs = std::abs(matrix_copy[i][j]);
                            max_i = i;
                            max_j = j;
                        }
                    }
                }
                break;
            }
            case model::EliminationChoiceMethod::Cyclic: {
                if (max_j == matrix_.size() - 1) {
                    max_i = (max_i + 1) % (matrix_.size() - 1);
                    max_j = max_i + 1;
                } else {
                    ++max_j;
                }
                break;
            }
            default:
                break;
        }

        non_diagonal_sum -= 2 * matrix_copy[max_i][max_j] * matrix_copy[max_i][max_j];

        double x = -2 * matrix_copy[max_i][max_j];
        double y = matrix_copy[max_i][max_i] - matrix_copy[max_j][max_j];
        double sine, cosine;
        if (y == 0) {
            cosine = 1 / std::sqrt(2);
            sine = cosine;
        } else {
            double root = std::sqrt(x * x + y * y);
            double value = std::abs(y) / root;
            cosine = std::sqrt(0.5 + 0.5 * value);
            sine = std::abs(x) / (2 * cosine * root);
            if (x * y < 0) {
                sine = -sine;
            }
        }

        matrix_copy.MultiplyByOrthogonalLeft(max_i, max_j, cosine, sine);
        matrix_copy.MultiplyByOrthogonalRight(max_i, max_j, cosine, -sine);
    } while (non_diagonal_sum >= eps);

    std::vector<double> eigenvalues(matrix_.size());
    for (std::size_t i = 0; i != matrix_copy.Size(); ++i) {
        eigenvalues[i] = matrix_copy[i][i];
    }

    return {std::move(eigenvalues), iteration_num};
}

bool Matrix::IsDiagonal() const {
    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        for (std::size_t j = 0; j != matrix_.size(); ++j) {
            if (matrix_[i][j] != 0 && i != j) {
                return false;
            }
        }
    }

    return true;
}

bool Matrix::IsSymmetric() const {
    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        for (std::size_t j = i + 1; j != matrix_.size(); ++j) {
            if (matrix_[i][j] != matrix_[j][i]) {
                return false;
            }
        }
    }

    return true;
}

bool Matrix::IsDiagonallyDominant() const {
    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        double sum = 0;
        for (std::size_t j = 0; j != matrix_.size(); ++j) {
            if (i != j) {
                sum += std::abs(matrix_[i][j]);
            }
        }
        if (std::abs(matrix_[i][i]) < sum) {
            return false;
        }
    }

    return true;
}

bool Matrix::IsTridiagonal() const {
    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        for (int j = 0; j != matrix_.size(); ++j) {
            if (matrix_[i][j] != 0 && i != j && i != j - 1 && i != j + 1) {
                return false;
            }
        }
    }

    return true;
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
    Matrix matrix(n);
    for (std::size_t i = 0; i != n; ++i) {
        matrix[i][i] = k;
    }
    return matrix;
}

Matrix Matrix::CreateGilbert(std::size_t n) {
    Matrix matrix(n);
    for (std::size_t i = 0; i != n; ++i) {
        for (std::size_t j = 0; j != n; ++j) {
            matrix[i][j] = 1 / (double)(i + j + 1);
        }
    }
    return matrix;
}

Matrix Matrix::CreateTridiagonal(std::size_t n) {
    Matrix matrix(n);
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

void Matrix::MakeSymmetric() {
    for (std::size_t i = 0; i != matrix_.size(); ++i) {
        for (std::size_t j = i + 1; j != matrix_.size(); ++j) {
            matrix_[i][j] = matrix_[j][i];
        }
    }
}

Matrix Matrix::CreateRandomDiagonallyDominant(bool is_symmetric, std::size_t n, double lower_bound,
                                              double upper_bound) {
    Matrix matrix(n);
    std::random_device r;
    std::default_random_engine re(r());
    std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    std::mt19937 rng(r());
    std::uniform_int_distribution<std::mt19937::result_type> dist(0, n * n - 1);

    for (std::size_t i = 0; i != n; ++i) {
        matrix[i][i] = unif(re);
    }

    for (std::size_t index = 0; index != n; ++index) {
        std::size_t i = 0;
        std::size_t j = 0;
        std::size_t num = 0;
        while (i == j || matrix[i][j] != 0) {
            num = dist(rng);
            i = num / n;
            j = num % n;
        }
        double value = 0;
        double sum = 0;
        for (std::size_t k = 0; k != n; ++k) {
            if (k != i) sum += std::abs(matrix[i][k]);
        }
        double gap = std::abs(matrix[i][i]) - sum;

        if (is_symmetric) {
            double second_sum = 0;
            for (std::size_t k = 0; k != n; ++k) {
                if (k != j) second_sum += std::abs(matrix[j][k]);
            }
            double second_gap = std::abs(matrix[j][j]) - second_sum;
            gap = std::min(gap, second_gap);
        }

        std::uniform_real_distribution<double> unif(-std::abs(gap) / 2, std::abs(gap) / 2);
        value = unif(re);
        matrix[i][j] = value;
        if (is_symmetric) {
            matrix[j][i] = value;
        }
    }

    return matrix;
}

Matrix Matrix::CreateRandomDiagonallyDominant(std::size_t n, double lower_bound,
                                              double upper_bound) {
    return CreateRandomDiagonallyDominant(false, n, lower_bound, upper_bound);
}

Matrix Matrix::CreateRandomSymmetric(std::size_t n, double lower_bound, double upper_bound) {
    Matrix matrix = CreateRandom(n, lower_bound, upper_bound);
    matrix.MakeSymmetric();
    return matrix;
}

Matrix Matrix::CreateRandomSymmetricDiagonallyDominant(std::size_t n, double lower_bound,
                                                       double upper_bound) {
    return CreateRandomDiagonallyDominant(true, n, lower_bound, upper_bound);
}

Matrix Matrix::CreateRandom(std::size_t n, double lower_bound, double upper_bound) {
    Matrix matrix(n);

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
