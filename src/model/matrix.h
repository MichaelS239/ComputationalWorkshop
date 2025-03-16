#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "util/matrix_util.h"

namespace model {
class Matrix {
private:
    std::vector<std::vector<double>> matrix_;

    void CalculateDeterminant(std::vector<bool> available_indices, unsigned depth,
                              double& result) const;

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

    std::size_t Size() const {
        return matrix_.size();
    }

    std::vector<std::vector<double>> GetMatrix() const {
        return matrix_;
    }

    double Determinant() const;
    double Norm() const;
    Matrix Inverse() const;
    std::vector<double> SolveSystem(std::vector<double> const& vector) const;
    double NormConditionNumber() const;
    double VolumeConditionNumber() const;
    double AngleConditionNumber() const;

    std::string ToString() const;

    static Matrix CreateDiagonal(std::size_t n, double k);
    static Matrix CreateGilbert(std::size_t n);
    static Matrix CreateTridiagonal(std::size_t n);
    static Matrix CreateRandom(std::size_t n, double lower_bound = -10, double upper_bound = 10);
};

}  // namespace model
