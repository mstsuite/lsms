
set(CMAKE_CXX_COMPILER_LAUNCHER ccache)
set(CMAKE_C_COMPILER_LAUNCHER ccache)

set(CMAKE_Fortran_FLAGS "-Og -g -fcheck=all -fbacktrace")
set(CMAKE_CXX_FLAGS "-Og -g")

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
set(BLA_VENDOR Generic)