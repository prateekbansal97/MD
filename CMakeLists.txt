#cmake_minimum_required(VERSION 3.18) # 3.18+ recommended for CUDA_ARCHITECTURES
#project(MD LANGUAGES CXX)
#
#set(CMAKE_CXX_STANDARD 20)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
#
#option(ENABLE_CUDA "Build with CUDA acceleration" OFF)
#
## ---- AmberTopology library ----
#add_library(amber_topology
#        src/util/IO.cpp
#        src/ChemicalEntity/atom.cpp
#        src/ChemicalEntity/element.cpp
#        src/AmberTopology/AmberTopology.cpp
#        src/ChemicalEntity/molecule.cpp
#        src/Bonded/harmonicUB.cpp
#        src/Bonded/HarmonicImproper.cpp
#        src/Bonded/HarmonicBond.cpp
#        src/Bonded/HarmonicAngle.cpp
#        src/Bonded/CosineDihedral.cpp
#        src/Bonded/cmapgroup.cpp
#        src/NonBonded/LennardJones.cpp
#        src/NonBonded/CoulombicEE.cpp
#        src/AmberTopology/PrintTopology.cpp
#        src/AmberTopology/CreateBondedFromParsed.cpp
#        src/AmberTopology/CreateNonBondedFromParsed.cpp
#        src/AmberTopology/TopologyParser.cpp
#)
#
#target_include_directories(amber_topology PUBLIC ${CMAKE_SOURCE_DIR}/include)
#target_compile_features(amber_topology PUBLIC cxx_std_17)
#
## ---- System library (your integrator/forces orchestration) ----
#add_library(md_system
#        src/System/System.cpp
#        src/System/metrics.cpp
#        include/System/metrics.h
#        src/System/init.cpp
#        src/System/NonBondedUtils.cpp
#        src/System/PBC.cpp
#        src/System/BondedEnergy.cpp
#        src/System/NonBondedEnergy.cpp
#        src/System/BondedForces.cpp
#        src/System/NonBondedForces.cpp
#)
#target_include_directories(md_system PUBLIC ${CMAKE_SOURCE_DIR}/include)
#target_link_libraries(md_system PUBLIC amber_topology)
#
## ---- Optional CUDA nonbonded library ----
#if(ENABLE_CUDA)
#    add_library(nonbonded_cuda
#            src/cuda/nonbonded_cuda.cu
#    )
#
#    target_include_directories(nonbonded_cuda PUBLIC ${CMAKE_SOURCE_DIR}/include)
#    target_compile_features(nonbonded_cuda PUBLIC cxx_std_17)
#
#    # Choose GPU architectures (edit these for your cluster)
#    set_target_properties(nonbonded_cuda PROPERTIES
#            CUDA_SEPARABLE_COMPILATION ON
#            CUDA_ARCHITECTURES 70;75;80;86
#    )
#
#    # Link CUDA runtime
#    find_package(CUDAToolkit REQUIRED)
#    target_link_libraries(nonbonded_cuda PUBLIC CUDA::cudart)
#
#    target_compile_definitions(md_system PUBLIC MD_ENABLE_CUDA=1)
#    target_link_libraries(md_system PUBLIC nonbonded_cuda)
#endif()
#
#
#find_package(PkgConfig REQUIRED)
#pkg_check_modules(FFTW3 REQUIRED fftw3)
#include_directories(${FFTW3_INCLUDE_DIRS})
#
## ---- Executable ----
#add_executable(main apps/main.cpp)
#target_link_libraries(main PRIVATE md_system)
#target_link_libraries(main ${FFTW3_LIBRARIES})

cmake_minimum_required(VERSION 3.18) # 3.18+ recommended for CUDA_ARCHITECTURES
project(MD LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ENABLE_CUDA "Build with CUDA acceleration" OFF)
# ---- Dependencies (FFTW3) ----
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)

# ---- AmberTopology library ----
add_library(amber_topology
        src/util/IO.cpp
        src/ChemicalEntity/atom.cpp
        src/ChemicalEntity/element.cpp
        src/AmberTopology/AmberTopology.cpp
        src/ChemicalEntity/molecule.cpp
        src/Bonded/harmonicUB.cpp
        src/Bonded/HarmonicImproper.cpp
        src/Bonded/HarmonicBond.cpp
        src/Bonded/HarmonicAngle.cpp
        src/Bonded/CosineDihedral.cpp
        src/Bonded/cmapgroup.cpp
        src/NonBonded/LennardJones.cpp
        src/NonBonded/CoulombicEE.cpp
        src/AmberTopology/PrintTopology.cpp
        src/AmberTopology/CreateBondedFromParsed.cpp
        src/AmberTopology/CreateNonBondedFromParsed.cpp
        src/AmberTopology/TopologyParser.cpp
)

target_include_directories(amber_topology PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_compile_features(amber_topology PUBLIC cxx_std_20)

# ---- System library (Integrator/Forces/PME) ----
add_library(md_system
        src/System/System.cpp
        src/System/metrics.cpp
        include/System/metrics.h
        src/System/init.cpp
        src/System/NonBondedUtils.cpp
        src/System/PBC.cpp
        src/System/BondedEnergy.cpp
        src/System/NonBondedEnergy.cpp
        src/System/BondedForces.cpp
        src/System/NonBondedForces.cpp
)
target_include_directories(md_system PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_link_libraries(md_system PUBLIC
        amber_topology
        PkgConfig::FFTW3
)

# ---- Optional CUDA nonbonded library ----
if(ENABLE_CUDA)
    add_library(nonbonded_cuda
            src/cuda/nonbonded_cuda.cu
    )

    target_include_directories(nonbonded_cuda PUBLIC ${CMAKE_SOURCE_DIR}/include)
    target_compile_features(nonbonded_cuda PUBLIC cxx_std_20)

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
add_executable(main apps/main.cpp)
# main only needs to link md_system; it inherits FFTW linkage from md_system
target_link_libraries(main PRIVATE md_system)