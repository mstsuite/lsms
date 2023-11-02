# Toolchain for building LSMS on MacOS

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)

set(MST_LINEAR_SOLVER_DEFAULT 0x0005)
set(MST_BUILD_KKR_MATRIX_DEFAULT 0x1000)
set(BLA_VENDOR Apple)

set(CMAKE_CXX_COMPILER "mpic++")
set(CMAKE_C_COMPILER "gcc-13")

set(BUILD_WITH_OPENMP OFF)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-O3 -mtune=native")
set(LINK_FLAGS "-Wl,-ld_classic")
# set(LINK_FLAGS "-ld64")
# set(libxc_LIBRARIES "/Users/vie/Documents/MuST/external/libxc/LibXC/lib/libxc.a")
# set(libxc_INCLUDE_DIR "/Users/vie/Documents/MuST/external/libxc/LibXC/include")
