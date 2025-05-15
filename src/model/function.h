#pragma once

#include <array>
#include <cstddef>
#include <functional>

namespace model {
using Func = std::function<double(double)>;
using TwoVariableFunc = std::function<double(double, double)>;

template <std::size_t N>
using MultivariableFunc = std::function<double(std::array<double, N> const&)>;
template <std::size_t N>
using GradientFunc = std::function<std::array<double, N>(std::array<double, N>)>;
template <std::size_t N>
using SecondDerivativeFunc =
        std::function<std::array<std::array<double, N>, N>(std::array<double, N>)>;
}  // namespace model
