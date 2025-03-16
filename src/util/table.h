#pragma once

#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

namespace util {
void PrintTable(std::vector<std::vector<double>> const& table,
                std::vector<std::string> const& headers);
void PrintTable(std::vector<std::pair<double, double>> const& table,
                std::pair<std::string, std::string> const& headers = {"x", "f(x)"});

template <typename T>
std::vector<std::pair<double, double>> CreateTable(T&& first, T&& second) {
    assert(first.size() == second.size());

    std::vector<std::pair<double, double>> table;
    for (std::size_t i = 0; i != first.size(); ++i) {
        table.emplace_back(first[i], second[i]);
    }
    return table;
}

}  // namespace util
