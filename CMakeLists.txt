cmake_minimum_required(VERSION 3.16)
project(PLAME_LP)

set(CMAKE_CXX_STANDARD 17)

add_executable(PLAME_LP main.cpp solvelp.h)
add_definitions(-DILOUSESTL -DIL_STD)
# ILOUSESTL: iostream.h: No such file or directory
# IL_STD: cerr was not declared in this scope; did you mean std::cerr

include_directories(
        /mnt/d/CPLEXLinux/concert/include
        /mnt/d/CPLEXLinux/cplex/include
) # WSL format

find_library(ILOCPLEX ilocplex /mnt/d/CPLEXLinux/cplex/lib/x86-64_linux/static_pic/)
find_library(CPLEX cplex /mnt/d/CPLEXLinux/cplex/lib/x86-64_linux/static_pic/)
find_library(CONCERT concert /mnt/d/CPLEXLinux/concert/lib/x86-64_linux/static_pic/)

target_link_libraries(PLAME_LP PUBLIC
        ${ILOCPLEX}
        ${CONCERT}
        ${CPLEX}
        libm.so # should use dynamic lib
        libpthread.so # should use dynamic lib
        ${CMAKE_DL_LIBS})

