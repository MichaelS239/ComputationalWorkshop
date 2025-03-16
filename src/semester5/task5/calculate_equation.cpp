#include "calculate_equation.h"

#include <cmath>
#include <utility>
#include <vector>

namespace semester5_task5 {

double Solution(double x) {
    return std::exp(-x) * (x + 1);
}

// y(0) = 1, y'(0) = 0, y''(0) = -1, y_3(0) = 2, y_4(0) = -3, y_5(0) = 4, y_6(0) = -5
double Taylor(double x) {
    return 1 - std::pow(x, 2) / 2 + std::pow(x, 3) / 3 - std::pow(x, 4) / 8 + std::pow(x, 5) / 30 -
           std::pow(x, 6) / 144;
}

double f(double x, double y) {
    return -y + std::exp(-x);
}

std::vector<double> CalculateSolution(int start, int end, double h) {
    std::vector<double> v;
    for (int i = start; i != end; ++i) {
        v.push_back(Solution(i * h));
    }
    return v;
}

std::vector<double> CalculateTaylor(int start, int end, double h) {
    std::vector<double> v;
    for (int i = start; i != end; ++i) {
        v.push_back(Taylor(i * h));
    }
    return v;
}

std::vector<double> CalculateAdams(int start, int end, double h) {
    start -= 5;
    std::vector<std::vector<double>> diffs =
            std::vector<std::vector<double>>(end - start, std::vector<double>(7));
    std::vector<double> table = CalculateTaylor(start, start + 5, h);
    for (int i = 0; i != table.size(); ++i) {
        diffs[i][0] = (start + i) * h;
        diffs[i][1] = table[i];
        diffs[i][2] = h * f(diffs[i][0], diffs[i][1]);
    }
    for (int i = 3; i != 7; ++i) {
        for (int j = 0; j != table.size() - i + 2; ++j) {
            diffs[j][i] = diffs[j + 1][i - 1] - diffs[j][i - 1];
        }
    }
    for (int i = table.size(); i != end - start; ++i) {
        diffs[i][0] = (start + i) * h;
        double sum = diffs[i - 1][2] + diffs[i - 2][3] / 2 + diffs[i - 3][4] * 5 / 12 +
                     diffs[i - 4][5] * 3 / 8 + diffs[i - 5][6] * 251 / 720;
        diffs[i][1] = diffs[i - 1][1] + sum;
        diffs[i][2] = h * f(diffs[i][0], diffs[i][1]);
        for (int j = 3; j != 7; ++j) {
            diffs[i - j + 2][j] = diffs[i - j + 3][j - 1] - diffs[i - j + 2][j - 1];
        }
    }

    std::vector<double> result;
    for (int i = start + 5; i != end; ++i) {
        result.push_back(diffs[i - start][1]);
    }
    return result;
}

std::vector<double> CalculateRungeKutta(int start, int end, double h) {
    std::vector<double> result;
    double x = 0;
    double y = 1;
    double k1, k2, k3, k4;
    double y1;
    for (int i = start; i != end; ++i) {
        k1 = h * f(x, y);
        k2 = h * f(x + h / 2, y + k1 / 2);
        k3 = h * f(x + h / 2, y + k2 / 2);
        k4 = h * f(x + h, y + k3);
        y1 = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        result.push_back(y1);
        x = i * h;
        y = y1;
    }
    return result;
}

std::vector<double> CalculateEuler(int start, int end, double h) {
    std::vector<double> result;
    double x = 0;
    double y = 1;
    double y1;
    for (int i = start; i != end; ++i) {
        y1 = y + h * f(x, y);
        result.push_back(y1);
        x = i * h;
        y = y1;
    }
    return result;
}

std::vector<double> CalculateEuler1(int start, int end, double h) {
    std::vector<double> result;
    double x = 0;
    double y = 1;
    double y1;
    for (int i = start; i != end; ++i) {
        y1 = y + h * f(x + h / 2, y + f(x, y) * h / 2);
        result.push_back(y1);
        x = i * h;
        y = y1;
    }
    return result;
}

std::vector<double> CalculateEuler2(int start, int end, double h) {
    std::vector<double> result;
    double x = 0;
    double y = 1;
    double y1;
    for (int i = start; i != end; ++i) {
        y1 = y + (f(x, y) + f(i * h, y + h * f(x, y))) * h / 2;
        result.push_back(y1);
        x = i * h;
        y = y1;
    }
    return result;
}

std::vector<std::vector<double>> CompareWithSolution(int start, int end, double h,
                                                     CalcMethod method) {
    std::vector<std::vector<double>> table;
    std::vector<double> result = method(start, end, h);
    for (int i = start; i != end; ++i) {
        std::vector<double> cur(4);
        cur[0] = i * h;
        cur[1] = Solution(i * h);
        cur[2] = result[i - start];
        cur[3] = std::abs(cur[1] - cur[2]);
        table.push_back(std::move(cur));
    }
    return table;
}

}  // namespace semester5_task5
