# Toolchain for building LSMS with HIP on OLCF Frontier with amd compiler

set(SEARCH_LAPACK OFF)
set(SEARCH_BLAS OFF)
set(BLAS_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_amd.a")
set(LAPACK_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_amd.a")
# set(USE_ACCELERATOR_HIP ON)
#set(MST_LINEAR_SOLVER_DEFAULT 0x0013)
set(MST_LINEAR_SOLVER_DEFAULT 0x0005)
#set(MST_BUILD_KKR_MATRIX_DEFAULT 0x3000)
set(MST_BUILD_KKR_MATRIX_DEFAULT 0x1000)

# Currently Loaded Modules:
#   1) craype-x86-trento                       6) cray-pmi/6.0.16      11) netlib-scalapack/2.1.0  16) cray-mpich/8.1.12
#   2) libfabric/1.13.0.0                      7) cray-pmi-lib/6.0.16  12) cmake/3.21.3            17) cray-libsci/21.08.1.2
#   3) craype-network-ofi                      8) DefApps/default      13) amd/4.5.0               18) PrgEnv-amd/8.2.0
#   4) perftools-base/21.12.0                  9) hdf5/1.12.0          14) craype/2.7.13
#   5) xpmem/2.3.1-2.2_26.4__g8e17d33.shasta  10) fftw/3.3.9           15) cray-dsmml/0.2.2
