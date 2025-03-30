#include "vector_util.h"

#include <cmath>
#include <cstddef>

namespace util {
double ScalarProduct(std::vector<double> const& vec1, std::vector<double> const& vec2) {
    double product = 0;
    for (std::size_t i = 0; i != vec1.size(); ++i) {
        product += vec1[i] * vec2[i];
    }
    return product;
}

double Norm(std::vector<double> const& vec) {
    return std::sqrt(ScalarProduct(vec, vec));
}

void Normalize(std::vector<double>& vec) {
    double norm = Norm(vec);
    for (std::size_t i = 0; i != vec.size(); ++i) {
        vec[i] /= norm;
    }
}

}  // namespace util
