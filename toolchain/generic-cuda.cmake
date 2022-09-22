#
# Toolchain for building LSMS on a generic Linux system with Cuda GPU
#

message(STATUS "Use toolchain file generic-cuda")

include_directories("/opt/amazon/openmpi/include/")

set(USE_ACCELERATOR_CUDA_C ON)

set(BUILD_TESTING OFF)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)

set(MST_LINEAR_SOLVER_DEFAULT 0x0013)
set(MST_BUILD_KKR_MATRIX_DEFAULT 0x3000)

set(CMAKE_CXX_COMPILER "mpic++")
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_Fortran_COMPILER "gfortran")

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-O3 -mtune=native")
set(CMAKE_Fortran_FLAGS "-O3 -mtune=native")
# set(CMAKE_OPTIMIZE_DEPENDENCIES TRUE)
# set(CMAKE_Fortran_PREPROCESS TRUE)
