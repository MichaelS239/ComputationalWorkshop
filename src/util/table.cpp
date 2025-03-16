#include "table.h"

#include <cstddef>
#include <iostream>
#include <sstream>
#include <vector>

namespace util {

void PrintTable(std::vector<std::vector<double>> const& table,
                std::vector<std::string> const& headers) {
    std::stringstream ss;
    ss.precision(15);
    std::vector<std::vector<std::size_t>> sizes;
    std::vector<std::size_t> max_sizes(headers.size(), 0);
    for (std::size_t i = 0; i != headers.size(); ++i) {
        max_sizes[i] = headers[i].size();
    }
    for (std::size_t i = 0; i != table.size(); ++i) {
        sizes.emplace_back();
        for (std::size_t j = 0; j != table[i].size(); ++j) {
            ss << table[i][j];
            std::size_t size = ss.str().size();
            ss.str("");
            ss.clear();
            max_sizes[j] = std::max(max_sizes[j], size);
            sizes[sizes.size() - 1].push_back(size);
        }
    }
    for (std::size_t i = 0; i != headers.size(); ++i) {
        std::cout << '+';
        for (std::size_t j = 0; j != max_sizes[i] + 2; ++j) std::cout << '-';
    }
    std::cout << "+\n";
    for (std::size_t i = 0; i != headers.size(); ++i) {
        std::cout << "| " << headers[i];
        for (std::size_t j = 0; j != max_sizes[i] - headers[i].size() + 1; ++j) std::cout << ' ';
    }
    std::cout << "|\n";
    for (std::size_t i = 0; i != headers.size(); ++i) {
        std::cout << '+';
        for (std::size_t j = 0; j != max_sizes[i] + 2; ++j) std::cout << '-';
    }
    std::cout << "+\n";
    for (std::size_t i = 0; i != table.size(); ++i) {
        for (std::size_t j = 0; j != table[i].size(); ++j) {
            ss << "| " << table[i][j];
            for (std::size_t k = sizes[i][j]; k != max_sizes[j] + 1; ++k) {
                ss << ' ';
            }
        }
        ss << "|\n";
        std::cout << ss.str();
        ss.str("");
        ss.clear();
    }

    for (std::size_t i = 0; i != headers.size(); ++i) {
        std::cout << '+';
        for (std::size_t j = 0; j != max_sizes[i] + 2; ++j) std::cout << '-';
    }
    std::cout << "+\n";
}

void PrintTable(std::vector<std::pair<double, double>> const& table,
                std::pair<std::string, std::string> const& headers) {
    std::vector<std::vector<double>> new_table(table.size());
    std::vector<std::string> new_headers = {headers.first, headers.second};
    for (std::size_t i = 0; i != table.size(); ++i) {
        new_table[i] = {table[i].first, table[i].second};
    }
    PrintTable(new_table, new_headers);
}

}  // namespace util
