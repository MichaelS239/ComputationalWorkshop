#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace util {
std::vector<double> InputPoints(std::size_t n);
char InputChoice(char first_choice, char second_choice);
int InputSegments(double a, double b, int minimum_number = 1);
double InputStep();
std::pair<double, double> InputBoundaries(std::string const& description);
}  // namespace util
