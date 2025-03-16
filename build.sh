#!/bin/bash

# Stop on error:
set -e

function print_help() {
cat << EOF
Usage: ./build.sh [options]

Possible options:
  -h,         --help                  Display help
  -s[N],      --semester[=N]          Set semester number N of a specific task to build (by default all tasks are built).
                                      Both semester and task number options must be set in order to build the specific task.
                                      Possible values of N: 5, 6.
  -t[N],      --task[=N]              Set task number N of a specific task to build (by default all tasks are built).
                                      Both semester and task number options must be set in order to build the specific task.
                                      Possible values of N: 1, 2, 3, 4.1, 4.2, 4.3, 5, 6 (if semester=5);
                                                            1 (if semester=6).
  -j[N],      --jobs[=N]              Allow N jobs at once (default [=1])
  -d,         --debug                 Set debug build type
EOF
}

for i in "$@"
    do
    case $i in
        -j*|--jobs=*) # Allow N jobs at once
            JOBS_OPTION=$i
            ;;
        -d|--debug) # Set debug build type
            DEBUG_MODE=true
            ;;
        --semester=*) # Set semester number, long option
            SEMESTER_NUM="${i#*=}"
            ;;
        -s*) # Set semester number, short option
            SEMESTER_NUM="${i#*s}"
            ;;
        --task=*) # Set task number, long option
            TASK_NUM="${i#*=}"
            ;;
        -t*) # Set task number, short option
            TASK_NUM="${i#*t}"
            ;;
        -h|--help|*) # Display help
            print_help
            exit 0
            ;;
    esac
done

mkdir -p lib
cd lib

if [[ ! -d "eigen" ]] ; then
  git clone https://gitlab.com/libeigen/eigen.git --branch 3.4.0 --depth 1
fi

if [[ $DEBUG_MODE == true ]]; then
  PREFIX="$PREFIX -D CMAKE_BUILD_TYPE=Debug"
fi

if [[ -n $SEMESTER_NUM ]]; then
  PREFIX="$PREFIX -D SEMESTER_NUM=${SEMESTER_NUM}"
fi

if [[ -n $TASK_NUM ]]; then
  PREFIX="$PREFIX -D TASK_NUM=${TASK_NUM}"
fi

cd ..
mkdir -p build
cd build
rm -f CMakeCache.txt
cmake $PREFIX .. && make $JOBS_OPTION
