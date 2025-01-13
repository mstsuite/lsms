# Toolchain for building LSMS with HIP on OLCF Frontier with amd compiler

set(CMAKE_BUILD_TYPE Debug)
set(BUILD_TESTING OFF)

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)

set(SEARCH_LAPACK OFF)
set(SEARCH_BLAS OFF)

set(BUILD_WITH_OPENMP OFF)

set(BLAS_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a")
set(LAPACK_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu.a")
set(USE_ACCELERATOR_HIP ON)
set(MST_LINEAR_SOLVER_DEFAULT 0x0020)
#set(MST_LINEAR_SOLVER_DEFAULT 0x0005)
set(MST_BUILD_KKR_MATRIX_DEFAULT 0x3000)
#set(MST_BUILD_KKR_MATRIX_DEFAULT 0x1000)
set(AMDGPU_TARGETS "gfx90a")
set(GPU_TARGETS "gfx90a")
# set(CMAKE_CXX_FLAGS "-std=c++17 -mtune=native")
# set(CMAKE_CXX_FLAGS "-std=c++17 --gcc-toolchain=/opt/cray/pe/gcc/12.2.0/snos -L/opt/cray/pe/gcc/12.2.0/snos/lib64 -mtune=native")

# Loaded Modules:
# module load PrgEnv-cray-amd
# module load cray-hdf5
# module load cmake libtool
# module load rocm

