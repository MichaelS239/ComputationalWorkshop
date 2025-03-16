#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>

#include "tasks.h"

void print_help() {
    std::cout << R"(Usage: ./ComputationalWorkshop [options]

Possible options:
-h,         --help                  Display help
-s[N],      --semester[=N]          Set semester number N of a specific task.
                                    Both semester and task number options must be set in order to run the specific task.
                                    Possible values of N: 5, 6.
-t[N],      --task[=N]              Set task number N of a specific task.
                                    Both semester and task number options must be set in order to run the specific task.
                                    Possible values of N: 1, 2, 3, 4.1, 4.2, 4.3, 5, 6 (if semester=5);
                                                          1 (if semester=6).)"
              << '\n';
}

int main(int argc, char** argv) {
    std::size_t semester_num = 0, task_num = 0;
    bool help = false;

    for (int i = 1; i != argc; ++i) {
        std::string opt(argv[i]);
        if (opt == "-h" || opt == "--help") {
            help = true;
            break;
        } else if (opt.substr(0, 2) == "-s" || opt.substr(0, 11) == "--semester=") {
            if (semester_num != 0) {
                throw std::runtime_error("Error: multiple declaration of semester number");
            } else {
                std::string number;
                if (opt.substr(0, 2) == "-s") {
                    number = opt.substr(2, opt.size() - 2);
                } else {
                    number = opt.substr(11, opt.size() - 11);
                }
                std::size_t pos;
                int num = std::stoi(number, &pos);
                if (pos == number.size() && (num == 5 || num == 6)) {
                    semester_num = num;
                } else {
                    throw std::runtime_error("Error: wrong semester number");
                }
            }
        } else if (opt.substr(0, 2) == "-t" || opt.substr(0, 7) == "--task=") {
            if (task_num != 0) {
                throw std::runtime_error("Error: multiple declaration of semester number");
            } else {
                std::string number;
                if (opt.substr(0, 2) == "-t") {
                    number = opt.substr(2, opt.size() - 2);
                } else {
                    number = opt.substr(7, opt.size() - 7);
                }
                if (number == "4.1") {
                    task_num = 7;
                } else if (number == "4.2") {
                    task_num = 8;
                } else if (number == "4.3") {
                    task_num = 9;
                } else {
                    std::size_t pos;
                    int num = std::stoi(number, &pos);
                    if (pos == number.size() && num >= 1 && num <= 6) {
                        task_num = num;
                    } else {
                        throw std::runtime_error("Error: wrong task number");
                    }
                }
            }
        }
    }

    if (help) {
        print_help();
        return 0;
    }

    if (semester_num == 0) {
        throw std::runtime_error("Error: semester number is not set");
    }
    if (task_num == 0) {
        throw std::runtime_error("Error: task number is not set");
    }
    semester_num -= 5;
    task_num -= 1;

    if ((semester_num == 0 && (task_num == 3 || task_num >= tasks::semester5_tasks.size())) ||
        (semester_num == 1 && task_num >= tasks::semester6_tasks.size())) {
        throw std::runtime_error("Error: wrong task number");
    }

    if (semester_num == 0) {
        tasks::semester5_tasks[task_num]();
    } else {
        tasks::semester6_tasks[task_num]();
    }

    return 0;
}
