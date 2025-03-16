# About
This repository contains tasks for the computational workshop. Each task shows the work of a numerical method or a group of methods.

The project covers the following numerical methods:
- Root-finding algorithms
- Polynomial interpolation using Lagrange polynomial
- Numerical differentiation
- Numerical integration
- Numerical methods for solving initial value problems for ODE

The topics of tasks are the following:
- Semester 5:
  - Task 1: Root-finding algorithms
  - Task 2: Lagrange interpolation
  - Task 3: Numerical differentiation
  - Task 4.1: Numerical weighted integral calculation
  - Task 4.2: Numerical integration using simple quadrature formulas
  - Task 4.3: Numerical integration using composite quadrature formulas
  - Task 5: Numerical methods for solving initial value problems for ODE
  - Task 6: Numerical weighted integral calculation using quadrature formulas with maximum order of accuracy
- Semester 6:
  - Task 1: Condition numbers of a matrix

# Prerequisites
The following software is required to build the project:
- GNU GCC, version 11+
- CMake, version 3.10+

# Build instructions
Firstly, clone the repository and change the current directory to the project directory.

In order to build the whole project, run the following command in the terminal:
```sh
./build.sh
```
In order to build the specific task, run the following command:
```sh
./build.sh --semester=$SEMESTER --task=$TASK
```
where `$SEMESTER` is the number of semester (5 or 6) and `$TASK` is the number of desired task.
For example, the following command builds the task 4.1 of semester 5:
```sh
./build.sh --semester=5 --task=4.1
```

# Usage
After the build, the executable file `./ComputationalWorkshop` is located in `path/to/project/build/target` directory.
Firstly, navigate to this directory:
```sh
cd build/target
```
If you have built the specific task, run the following command to execute the task:
```sh
./ComputationalWorkshop
```
If you have built the whole project, you need to specify the task to run:
```sh
./ComputationalWorkshop --semester=$SEMESTER --task=$TASK
```
where `$SEMESTER` is the number of semester and `$TASK` is the number of desired task.
