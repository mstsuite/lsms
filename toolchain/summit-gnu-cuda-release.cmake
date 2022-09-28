#
# Toolchain for building LSMS with CUDA on OLCF Summit
#
# module load gcc/9.3.0  spectrum-mpi/10.4.0.3-20210112 hdf5/1.12.1 essl cuda ninja cmake
#

message(STATUS "Use toolchain file SUMMIT")

set(BUILD_TESTING ON)

set(USE_ACCELERATOR_CUDA_C ON)
set(USE_ESSL ON)
set(ARCH_IBM ON)
set(MST_LINEAR_SOLVER_DEFAULT 0x0013)
set(MST_BUILD_KKR_MATRIX_DEFAULT 0x3000)

set(CMAKE_CXX_COMPILER "mpic++")
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_Fortran_COMPILER "gfortran")

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-O3 -mtune=native -mcpu=native")
set(CMAKE_Fortran_FLAGS "-O3 -mtune=native -mcpu=native")
set(CMAKE_Fortran_PREPROCESS TRUE)



