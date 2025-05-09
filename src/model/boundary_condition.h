#pragma once

namespace model {

enum class BoundaryConditionKind { FirstKind, SecondKind, ThirdKind };

/**
 * Boundary condition of first kind:
 * u(a) = A
 * u(b) = B
 * Boundary condition of second kind:
 * u'(a) = A
 * u'(b) = B
 */
template <BoundaryConditionKind T>
struct BoundaryCondition {
    // a
    double left_boundary;
    // b
    double right_boundary;
    // A
    double left_value;
    // B
    double right_value;
};

/**
 * Boundary condition of third kind:
 * a0 * u(a) + a1 * u'(b) = A
 * b0 * u(b) + b1 * u'(b) = B
 */
template <>
struct BoundaryCondition<BoundaryConditionKind::ThirdKind> {
    // a
    double left_boundary;
    // a0
    double left_value_coef;
    // a1
    double left_derivative_coef;
    // A
    double left_value;
    // b
    double right_boundary;
    // b0
    double right_value_coef;
    // b1
    double right_derivative_coef;
    // B
    double right_value;
};

}  // namespace model
