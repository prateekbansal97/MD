cmake_minimum_required(VERSION 3.10)
project(MD)

set(CMAKE_C_COMPILER "/usr/bin/clang")
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(HEADERS io.h atom.h topology.h element.h harmonicUB.h HarmonicImproper.h HarmonicBond.h HarmonicAngle.h CosineDihedral.h CMapGroup.h)
set(SOURCES main.cpp io.cpp atom.cpp element.cpp topology.cpp harmonicUB.cpp HarmonicImproper.cpp HarmonicBond.cpp HarmonicAngle.cpp CosineDihedral.cpp molecule.cpp)
add_executable(main ${SOURCES})
