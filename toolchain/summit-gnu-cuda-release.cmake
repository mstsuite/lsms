#
# Toolchain for building LSMS with CUDA on OLCF Summit
#
# Currently Loaded Modules:
#  1) gcc/9.1.0   2) nsight-compute/2021.2.1   3) nsight-systems/2021.3.1.54   
#  4) cuda/11.0.3   5) essl/6.3.0   6) spectrum-mpi/10.4.0.3-20210112   
#  7) cmake/3.23.1   8) hdf5/1.12.1
#

message(STATUS "Use toolchain file")

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
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
set(CMAKE_OPTIMIZE_DEPENDENCIES TRUE)
set(CMAKE_Fortran_PREPROCESS TRUE)


