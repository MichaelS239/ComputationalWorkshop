#pragma once

#include <vector>

namespace util {
double ScalarProduct(std::vector<double> const& vec1, std::vector<double> const& vec2);
double Norm(std::vector<double> const& vec);
void Normalize(std::vector<double>& vec);
}  // namespace util
