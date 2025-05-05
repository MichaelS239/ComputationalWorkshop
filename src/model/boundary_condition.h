#pragma once

namespace model {

enum class BoundaryConditionKind { FirstKind, SecondKind, ThirdKind };

template <BoundaryConditionKind T>
struct BoundaryCondition {
    double left_boundary;
    double right_boundary;
    double left_value;
    double right_value;
};

template <>
struct BoundaryCondition<BoundaryConditionKind::ThirdKind> {
    double left_boundary;
    double left_value_coef;
    double left_derivative_coef;
    double left_value;
    double right_boundary;
    double right_value_coef;
    double right_derivative_coef;
    double right_value;
};

}  // namespace model
