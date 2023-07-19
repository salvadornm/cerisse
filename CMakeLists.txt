cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(CNS 
        DESCRIPTION "An AMR Compressible Reactive Flow Simulation Code for CPU and GPU"
        LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/EB_CNS/Tools/CMake")

############################### OPTIONS ###############################
# Base options (Must be included in case CMakeLists)
set(CNS_DIM             "3"        CACHE STRING "Number of physical dimensions")
set(CNS_PRECISION       "DOUBLE"   CACHE STRING "Floating point precision SINGLE or DOUBLE")
# set(CNS_CHEMISTRY_MODEL "LiDryer"  CACHE STRING "Chemical mechanism")
# set(CNS_EOS_MODEL       "Fuego"    CACHE STRING "Equation of state model")
# set(CNS_TRANSPORT_MODEL "Simple"   CACHE STRING "Transport model")

# General Options 
option(CNS_ENABLE_EB           "Enable Embedded Boundary"   OFF)
option(CNS_ENABLE_FCOMPARE     "Enable building fcompare"   OFF)
option(CNS_ENABLE_HYPRE        "Enable Offload to Hypre"    OFF)
option(CNS_ENABLE_TINY_PROFILE "Enable AMReX tiny profiler" OFF)
option(CNS_ENABLE_HDF5         "Enable HDF5 PlotFiles"      OFF)
option(CNS_ENABLE_TESTING      "Enable testing using CTest" ON)

# HPC Options
option(CNS_ENABLE_MPI     "Enable MPI"    OFF)
option(CNS_ENABLE_OPENMP  "Enable OpenMP" OFF)
option(CNS_ENABLE_CUDA    "Enable CUDA"   OFF)
option(CNS_ENABLE_HIP     "Enable HIP"    OFF)
option(CNS_ENABLE_SYCL    "Enable SyCL"   OFF)

# C++ Options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Setup Languages (is done in amrex?)
if(CNS_ENABLE_MPI)
  find_package(MPI REQUIRED CXX)
endif()
if (CNS_ENABLE_CUDA)
  enable_language(CUDA)
  if(CMAKE_VERSION VERSION_LESS 3.20)
    include(AMReX_SetupCUDA)
  endif()
endif()
if (CNS_ENABLE_HIP)
  enable_language(HIP)
endif()

############################## AMReX ###############################

message(STATUS "AMReX Configure Section")
set(AMREX_DIR "${CMAKE_SOURCE_DIR}/Submodules/amrex")
include(SetAMReXOptions) # This must be before add_subdirectory
list(APPEND CMAKE_MODULE_PATH "${AMREX_DIR}/Tools/CMake")
add_subdirectory(${AMREX_DIR})

############################# SUNDIALS #############################

message(STATUS "SUNDIALS Configure Section")
set(SUNDIALS_DIR "${CMAKE_SOURCE_DIR}/Submodules/sundials")
include(SetSundialsOptions)
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
add_subdirectory(${SUNDIALS_DIR})

########################### Print Summary ##########################

message(STATUS "CMAKE Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

# include(SetRpath)

add_subdirectory(EB_CNS/Exec) # This must be after all AMReX CMakes

# Testing
if(CNS_ENABLE_TESTING)
  enable_testing()
  include(CTest)
  add_subdirectory(EB_CNS/Tools/CTest)
endif()