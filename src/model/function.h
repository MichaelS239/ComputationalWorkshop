#pragma once

#include <functional>

namespace model {
using Func = std::function<double(double)>;
using TwoVariableFunc = std::function<double(double, double)>;
}  // namespace model
