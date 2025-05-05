#pragma once

namespace model {

enum class BoundaryConditionKind { FirstKind, SecondKind };

struct BoundaryCondition {
    double left_boundary;
    double right_boundary;
    double left_value;
    double right_value;
    BoundaryConditionKind kind;
};

}  // namespace model
