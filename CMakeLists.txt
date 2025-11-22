cmake_minimum_required(VERSION 3.10)
project(MD)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(HEADERS io.h atom.h)
set(SOURCES main.cpp io.cpp atom.cpp)
add_executable(main ${SOURCES})