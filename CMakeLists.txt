cmake_minimum_required(VERSION 3.10)
project(MD)

set(CMAKE_C_COMPILER "/usr/bin/clang")
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
#include_directories(${CMAKE_SOURCE_DIR}/include/AmberTopology/)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(HEADERS include/AmberTopology/io.h include/AmberTopology/atom.h include/AmberTopology/topology.h include/AmberTopology/element.h include/AmberTopology/harmonicUB.h include/AmberTopology/HarmonicImproper.h include/AmberTopology/HarmonicBond.h include/AmberTopology/HarmonicAngle.h include/AmberTopology/CosineDihedral.h include/AmberTopology/CMapGroup.h)
set(SOURCES main.cpp src/AmberTopology/io.cpp src/AmberTopology/atom.cpp src/AmberTopology/element.cpp src/AmberTopology/topology.cpp src/AmberTopology/harmonicUB.cpp src/AmberTopology/HarmonicImproper.cpp src/AmberTopology/HarmonicBond.cpp src/AmberTopology/HarmonicAngle.cpp src/AmberTopology/CosineDihedral.cpp src/AmberTopology/molecule.cpp
        src/System/system.cpp
        include/System/System.h)
add_executable(main ${SOURCES})
