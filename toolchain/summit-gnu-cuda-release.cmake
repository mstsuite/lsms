#
# Toolchain for building LSMS with CUDA on OLCF Summit
#
# module load gcc/9.3.0  spectrum-mpi/10.4.0.3-20210112 hdf5/1.12.1 essl cuda ninja cmake
#

message(STATUS "Use toolchain file SUMMIT")

set(BUILD_TESTING OFF)

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

# Currently Loaded Modules:
#   1) lsf-tools/2.0       6) cmake/3.27.7                    11) nsight-compute/2023.2.2
#   2) hsi/5.0.2.p5        7) spectrum-mpi/10.4.0.6-20230210  12) cuda/11.7.1
#   3) xalt/1.2.1          8) fftw/3.3.10                     13) hdf5/1.14.3
#   4) DefApps             9) netlib-scalapack/2.2.0
#   5) essl/6.1.0-1  (H)  10) gcc/9.3.0-compiler_only



