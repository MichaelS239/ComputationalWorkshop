#include "polynomial.h"

#include <cstddef>
#include <random>
#include <vector>

namespace model {
Polynomial Polynomial::CreateFromRandom(std::size_t deg, double lower_bound, double upper_bound) {
    Polynomial polynomial = std::vector<double>(deg, 0);
    std::random_device r;
    std::default_random_engine re(r());
    std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    for (std::size_t i = 0; i != deg; ++i) {
        polynomial[i] = unif(re);
    }

    return polynomial;
}

}  // namespace model
