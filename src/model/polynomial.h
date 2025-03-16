#pragma once

#include <cmath>
#include <cstddef>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace model {
class Polynomial {
private:
    std::vector<double> polynomial_;

public:
    Polynomial() = default;

    Polynomial(std::vector<double> const& polynomial) : polynomial_(polynomial) {}

    Polynomial(std::vector<double>&& polynomial) : polynomial_(std::move(polynomial)) {}

    std::vector<double>::size_type Degree() const {
        return polynomial_.size();
    }

    double operator()(double x) const {
        double s = 0;
        for (std::size_t i = 0; i != polynomial_.size(); ++i) {
            s += polynomial_[i] * std::pow(x, i);
        }
        return s;
    }

    double Prime(double x) const {
        double s = 0;
        for (std::size_t i = 1; i != polynomial_.size(); ++i) {
            s += polynomial_[i] * std::pow(x, i - 1) * i;
        }
        return s;
    }

    double Integral(double a, double b) const {
        double integral = 0;
        for (std::size_t i = 0; i != polynomial_.size(); ++i) {
            integral += polynomial_[i] / (i + 1) * (std::pow(b, i + 1) - std::pow(a, i + 1));
        }
        return integral;
    }

    double WeightedIntegral(std::vector<double> const& moments) const {
        double polynomial_integral = 0;
        for (std::size_t i = 0; i != polynomial_.size(); ++i) {
            polynomial_integral += polynomial_[i] * moments[i];
        }
        return polynomial_integral;
    }

    double& operator[](std::vector<double>::size_type i_e) {
        return polynomial_[i_e];
    }

    double const& operator[](std::vector<double>::size_type i_e) const {
        return polynomial_[i_e];
    }

    std::vector<double>::iterator begin() {
        return polynomial_.begin();
    }

    std::vector<double>::const_iterator begin() const {
        return polynomial_.begin();
    }

    std::vector<double>::iterator end() {
        return polynomial_.end();
    }

    std::vector<double>::const_iterator end() const {
        return polynomial_.end();
    }

    std::string ToString() const {
        std::stringstream ss;
        ss << "P(x) = ";
        for (std::size_t i = 0; i != polynomial_.size(); ++i) {
            if (i != 0) {
                if (polynomial_[i] < 0)
                    ss << " - ";
                else
                    ss << " + ";
            } else if (polynomial_[i] < 0)
                ss << "-";
            ss << std::abs(polynomial_[i]);
            if (i == 1)
                ss << " * x";
            else if (i != 0)
                ss << " * x^" << i;
        }
        return ss.str();
    }

    static Polynomial CreateFromRandom(std::size_t deg, double lower_bound = -10,
                                       double upper_bound = 10);
};

}  // namespace model
