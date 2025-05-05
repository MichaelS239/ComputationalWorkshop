#pragma once

#include <functional>
#include <vector>

#include "semester5/tasks.h"
#include "semester6/tasks.h"

namespace tasks {
std::vector<std::function<void(void)>> const semester5_tasks = {
        Semester5Task1, Semester5Task2,   Semester5Task3,   []() {},         Semester5Task5,
        Semester5Task6, Semester5Task4_1, Semester5Task4_2, Semester5Task4_3};
std::vector<std::function<void(void)>> const semester6_tasks = {Semester6Task1, Semester6Task2,
                                                                Semester6Task3, Semester6Task4,
                                                                Semester6Task5, Semester6Task6};

}  // namespace tasks
