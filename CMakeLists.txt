cmake_minimum_required(VERSION 3.10)

project(ComputationalWorkshop)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  "${CMAKE_BINARY_DIR}/target")

set(SEMESTER_NUM "-1" CACHE STRING "Semester number [5|6]")
set(TASK_NUM "-1" CACHE STRING "Task number")

include_directories(
    "src"
    "src/semester5"
    "src/semester6"
    "src/util"
    "src/model"
    "src/semester5/task1"
    "src/semester5/task2"
    "src/semester5/task3"
    "src/semester5/task4.1"
    "src/semester5/task4.2"
    "src/semester5/task4.3"
    "src/semester5/task5"
    "src/semester5/task6"
    "src/semester6/task1"
)

include_directories(SYSTEM "lib/eigen")

set(BINARY ${CMAKE_PROJECT_NAME})

if(SEMESTER_NUM STREQUAL "5")
    if(TASK_NUM STREQUAL "1" OR TASK_NUM STREQUAL "2" OR TASK_NUM STREQUAL "3" OR TASK_NUM STREQUAL "4.1" OR TASK_NUM STREQUAL "4.2" 
    OR TASK_NUM STREQUAL "4.3" OR TASK_NUM STREQUAL "5" OR TASK_NUM STREQUAL "6")
        file(GLOB_RECURSE lib_sources "src/util/*.h*" "src/util/*.cpp*" "src/util/*.cc*" "src/model/*.h*" "src/model/*.cpp*" "src/model/*.cc*" 
        src/semester5/task${TASK_NUM}/*.h* src/semester5/task${TASK_NUM}/*.cpp* src/semester5/task${TASK_NUM}/*.cc*)
        set(run_sources src/semester5/task${TASK_NUM}/main.cpp)
    elseif(TASK_NUM STREQUAL "-1")
        message(FATAL_ERROR "Error: no task number")
    else()
        message(FATAL_ERROR "Error: wrong task number")
    endif()
elseif(SEMESTER_NUM STREQUAL "6")
    if(TASK_NUM STREQUAL "1")
        file(GLOB_RECURSE lib_sources "src/util/*.h*" "src/util/*.cpp*" "src/util/*.cc*" "src/model/*.h*" "src/model/*.cpp*" "src/model/*.cc*"
        "src/semester6/task1/*.h*" "src/semester6/task1/*.cpp*" "src/semester6/task1/*.cc*")
        set(run_sources "src/semester6/task1/main.cpp")
    elseif(TASK_NUM STREQUAL "-1")
        message(FATAL_ERROR "Error: no task number")
    else()
        message(FATAL_ERROR "Error: wrong task number")
    endif()
elseif(SEMESTER_NUM STREQUAL "-1" AND TASK_NUM STREQUAL "-1")
    file(GLOB_RECURSE lib_sources "src/*.h*" "src/*.cpp*" "src/*.cc*")
    set(run_sources "src/main.cpp")
else()
    message(FATAL_ERROR "Error: no semester number")
endif()

add_library(${BINARY}_lib STATIC ${lib_sources})

add_executable(${BINARY} ${run_sources})
target_link_libraries(${BINARY} PRIVATE ${CMAKE_PROJECT_NAME}_lib)
