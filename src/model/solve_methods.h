#pragma once

namespace model {
enum class SolveMethod { Library, GaussElimination, LUDecomposition, QRDecomposition };

enum class CalculationMethod { Library, LUDecomposition, QRDecomposition };

enum class IterationMethod { Jacobi, Seidel };

enum class EigenvalueMethod { Power, ScalarProduct };

}  // namespace model
