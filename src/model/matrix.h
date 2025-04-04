#pragma once

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "model/solve_methods.h"

namespace model {
class Matrix {
private:
    class MatrixRow {
    private:
        double* row_;

    public:
        MatrixRow(double* row) : row_(row) {}

        double const& operator[](std::size_t i) const {
            return row_[i];
        }

        double& operator[](std::size_t i) {
            return row_[i];
        }
    };

    struct EigenInfo {
        double eigenvalue;
        std::vector<double> eigenvector;
        std::size_t iteration_num;
    };

    struct EigenValuesInfo {
        std::vector<double> eigenvalues;
        std::size_t iteration_num;
    };

    std::vector<std::vector<double>> matrix_;
    std::shared_ptr<Matrix> l_cache_;
    std::shared_ptr<Matrix> u_cache_;
    std::shared_ptr<Matrix> q_cache_;
    std::shared_ptr<Matrix> r_cache_;

    void CalculateDeterminant(std::vector<bool> available_indices, unsigned depth,
                              double& result) const;
    std::vector<double> SolveUpperTriangularSystem(std::vector<double> const& vector) const;
    std::vector<double> SolveLowerTriangularSystem(std::vector<double> const& vector) const;
    void MultiplyByOrthogonalLeft(std::size_t i, std::size_t j, double cosine, double sine);
    void MultiplyByOrthogonalRight(std::size_t i, std::size_t j, double cosine, double sine);

    std::pair<std::vector<double>, std::size_t> JacobiIteration(std::vector<double> const& vector,
                                                                double eps) const;
    std::pair<std::vector<double>, std::size_t> SeidelIteration(std::vector<double> const& vector,
                                                                double eps) const;
    void MakeSymmetric();
    static Matrix CreateRandomDiagonallyDominant(bool is_symmetric, std::size_t n,
                                                 double lower_bound, double upper_bound);
    EigenInfo PowerMethod(double eps) const;
    EigenInfo ScalarProductMethod(double eps) const;

    double CalculateNonDiagonalSum() const;

public:
    Matrix() = default;

    Matrix(std::vector<std::vector<double>> const& matrix) {
        for (std::size_t i = 0; i != matrix.size(); ++i) {
            if (matrix[i].size() != matrix.size()) {
                throw std::invalid_argument("Error: matrix is not square");
            }
        }

        matrix_ = matrix;
    }

    Matrix(std::vector<std::vector<double>>&& matrix) {
        for (std::size_t i = 0; i != matrix.size(); ++i) {
            if (matrix[i].size() != matrix.size()) {
                throw std::invalid_argument("Error: matrix is not square");
            }
        }

        matrix_ = std::move(matrix);
    }

    Matrix(std::size_t n) : matrix_(n, std::vector<double>(n)) {}

    Matrix(std::initializer_list<std::vector<double>> matrix) : matrix_(matrix) {}

    std::size_t Size() const {
        return matrix_.size();
    }

    std::vector<std::vector<double>> GetMatrix() const {
        return matrix_;
    }

    MatrixRow operator[](std::size_t i) {
        return &matrix_[i][0];
    }

    double Determinant(CalculationMethod calc_method = CalculationMethod::Library) const;
    double Norm() const;
    Matrix Inverse(SolveMethod calc_method = SolveMethod::Library) const;
    Matrix Transpose() const;
    void SelfTranspose();
    std::vector<double> MultiplyByVector(std::vector<double> const& vector) const;
    std::vector<double> SolveSystem(std::vector<double> const& vector,
                                    SolveMethod const solve_method = SolveMethod::Library) const;
    std::pair<std::vector<double>, std::size_t> SolveSystem(
            std::vector<double> const& vector, double eps,
            IterationMethod const solve_method = IterationMethod::Jacobi) const;
    std::pair<Matrix, std::vector<double>> GaussElimination(
            std::vector<double> const& vector) const;
    std::pair<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>> LUDecomposition() const;
    std::pair<std::shared_ptr<Matrix>, std::shared_ptr<Matrix>> QRDecomposition() const;
    void CacheLUDecomposition();
    void CacheQRDecomposition();

    double NormConditionNumber() const;
    double VolumeConditionNumber() const;
    double AngleConditionNumber() const;

    EigenInfo MaxAbsoluteEigenvalue(
            double eps, EigenvalueMethod const eigenvalue_method = EigenvalueMethod::Power) const;
    Matrix::EigenInfo MinAbsoluteEigenvalue(
            double eps, EigenvalueMethod const eigenvalue_method = EigenvalueMethod::Power) const;

    Matrix::EigenValuesInfo GetEigenValues(double eps) const;

    bool IsDiagonal() const;
    bool IsSymmetric() const;
    bool IsDiagonallyDominant() const;

    std::string ToString() const;

    static Matrix CreateDiagonal(std::size_t n, double k);
    static Matrix CreateGilbert(std::size_t n);
    static Matrix CreateTridiagonal(std::size_t n);
    static Matrix CreateRandomDiagonallyDominant(std::size_t n, double lower_bound = -10,
                                                 double upper_bound = 10);
    static Matrix CreateRandomSymmetric(std::size_t n, double lower_bound = -10,
                                        double upper_bound = 10);
    static Matrix CreateRandomSymmetricDiagonallyDominant(std::size_t n, double lower_bound = -10,
                                                          double upper_bound = 10);
    static Matrix CreateRandom(std::size_t n, double lower_bound = -10, double upper_bound = 10);
};

}  // namespace model
