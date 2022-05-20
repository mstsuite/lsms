# Toolchain for building LSMS with HIP on OLCF Frontier with amd compiler

set(USE_LIBXC OFF)
set(SEARCH_LAPACK OFF)
set(SEARCH_BLAS OFF)
set(BLAS_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu_82.a")
set(LAPACK_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_gnu_82.a")
# set(USE_ACCELERATOR_HIP ON)
#set(MST_LINEAR_SOLVER_DEFAULT 0x0013)
set(MST_LINEAR_SOLVER_DEFAULT 0x0001)
#set(MST_BUILD_KKR_MATRIX_DEFAULT 0x3000)
set(MST_BUILD_KKR_MATRIX_DEFAULT 0x1000)

# Currently Loaded Modules:
#   1) craype-x86-trento                       7) cray-pmi-lib/6.0.16  13) cray-libsci/21.08.1.2
#   2) libfabric/1.13.0.0                      8) cmake/3.22.1         14) PrgEnv-gnu/8.2.0
#   3) craype-network-ofi                      9) gcc/11.2.0           15) DefApps/default
#   4) perftools-base/21.12.0                 10) craype/2.7.13        16) hdf5/1.12.0
#   5) xpmem/2.3.2-2.2_1.16__g9ea452c.shasta  11) cray-dsmml/0.2.2     17) fftw/3.3.9
#   6) cray-pmi/6.0.16                        12) cray-mpich/8.1.12    18) netlib-scalapack/2.1.0
