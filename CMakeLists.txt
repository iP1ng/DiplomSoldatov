cmake_minimum_required(VERSION 3.7)
project(DiplomSoldatov)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp constants.h structures.h Gauss/gauss.cpp Gauss/gauss.h Triangle/triangle.cpp Triangle/triangle.h functions.h)
add_executable(DiplomSoldatov ${SOURCE_FILES})