#cmake_minimum_required(VERSION 3.10)
#project(MD)
#
#set(CMAKE_C_COMPILER "/usr/bin/clang")
#set(CMAKE_CXX_COMPILER "/usr/bin/clang++")
##include_directories(${CMAKE_SOURCE_DIR}/include/AmberTopology/)
#set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#
#set(HEADERS include/AmberTopology/io.h include/AmberTopology/atom.h include/AmberTopology/topology.h include/AmberTopology/element.h include/AmberTopology/harmonicUB.h include/AmberTopology/HarmonicImproper.h include/AmberTopology/HarmonicBond.h include/AmberTopology/HarmonicAngle.h include/AmberTopology/CosineDihedral.h include/AmberTopology/CMapGroup.h)
#set(SOURCES main.cpp src/AmberTopology/io.cpp src/AmberTopology/atom.cpp src/AmberTopology/element.cpp src/AmberTopology/topology.cpp src/AmberTopology/harmonicUB.cpp src/AmberTopology/HarmonicImproper.cpp src/AmberTopology/HarmonicBond.cpp src/AmberTopology/HarmonicAngle.cpp src/AmberTopology/CosineDihedral.cpp src/AmberTopology/molecule.cpp
#        src/System/system.cpp
#        include/System/System.h
#        src/AmberTopology/CMapGroup.cpp)
#add_executable(main ${SOURCES})


cmake_minimum_required(VERSION 3.18) # 3.18+ recommended for CUDA_ARCHITECTURES
project(MD LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ENABLE_CUDA "Build with CUDA acceleration" OFF)

# ---- AmberTopology library ----
add_library(amber_topology
        src/AmberTopology/io.cpp
        src/AmberTopology/atom.cpp
        src/AmberTopology/element.cpp
        src/AmberTopology/topology.cpp
        src/AmberTopology/harmonicUB.cpp
        src/AmberTopology/HarmonicImproper.cpp
        src/AmberTopology/HarmonicBond.cpp
        src/AmberTopology/HarmonicAngle.cpp
        src/AmberTopology/CosineDihedral.cpp
        src/AmberTopology/CMapGroup.cpp
        src/AmberTopology/molecule.cpp
        src/AmberTopology/LennardJones.cpp
        include/AmberTopology/LennardJones.h
        include/AmberTopology/CoulombicEE.h
        src/AmberTopology/CoulombicEE.cpp
)

target_include_directories(amber_topology PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_compile_features(amber_topology PUBLIC cxx_std_17)

# ---- System library (your integrator/forces orchestration) ----
add_library(md_system
        src/System/system.cpp
        src/System/metrics.cpp
        include/System/Metrics.h
)
target_include_directories(md_system PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(md_system PUBLIC amber_topology)

# ---- Optional CUDA nonbonded library ----
if(ENABLE_CUDA)
    add_library(nonbonded_cuda
            src/cuda/nonbonded_cuda.cu
    )

    target_include_directories(nonbonded_cuda PUBLIC ${CMAKE_SOURCE_DIR}/include)
    target_compile_features(nonbonded_cuda PUBLIC cxx_std_17)

    # Choose GPU architectures (edit these for your cluster)
    set_target_properties(nonbonded_cuda PROPERTIES
            CUDA_SEPARABLE_COMPILATION ON
            CUDA_ARCHITECTURES 70;75;80;86
    )

    # Link CUDA runtime
    find_package(CUDAToolkit REQUIRED)
    target_link_libraries(nonbonded_cuda PUBLIC CUDA::cudart)

    target_compile_definitions(md_system PUBLIC MD_ENABLE_CUDA=1)
    target_link_libraries(md_system PUBLIC nonbonded_cuda)
endif()

# ---- Executable ----
add_executable(main main.cpp)
target_link_libraries(main PRIVATE md_system)
