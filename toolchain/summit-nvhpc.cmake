# Toolchain for building LSMS with CUDA on OLCF Summit

# /sw/summit/nvhpc_sdk/Linux_ppc64le/23.9/compilers/include

set(ENV{INCLUDE} "$ENV{OLCF_NVHPC_ROOT}/compilers/include:$ENV{INCLUDE}")

set(USE_ACCELERATOR_CUDA_C ON)
set(USE_ESSL ON)
set(ARCH_IBM ON)
set(MST_LINEAR_SOLVER_DEFAULT 0x0013)
set(MST_BUILD_KKR_MATRIX_DEFAULT 0x3000)

