
message(STATUS "Use toolchain file for testing")

set(CMAKE_CXX_COMPILER_LAUNCHER ccache)
set(CMAKE_C_COMPILER_LAUNCHER ccache)

set(BUILD_TESTING ON)

set(MST_LINEAR_SOLVER_DEFAULT 0x0005)
set(MST_BUILD_KKR_MATRIX_DEFAULT 0x1000)

set(CMAKE_Fortran_FLAGS "-O2 -g -fcheck=bounds,mem,do -fbacktrace -Wall")
set(CMAKE_C_FLAGS "-O2 -g -Wall")
set(CMAKE_CXX_FLAGS "-O2 -g -Wall")

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
set(BLA_VENDOR Generic)
