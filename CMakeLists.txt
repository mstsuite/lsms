#
# Created by Franco P. Moitzi
#

cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

cmake_host_system_information(RESULT Ncpu QUERY NUMBER_OF_PHYSICAL_CORES)
message(STATUS "CMake ${CMAKE_VERSION} with max. ${Ncpu} threads")

project(lsms LANGUAGES Fortran C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# Default build type is RELEASE
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RELEASE CACHE BOOL
            "Choose the type of build, options are: Debug Release"
            FORCE)
endif (NOT CMAKE_BUILD_TYPE)

message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

# Build external libraries
message(STATUS "Build external libraries")

# Build LUA in this project
include(lua)

# Build LibXC in this project
include(libxc)

# Build mjson library
include(mjson)

#
# Enable Testing with CTest
#
option(BUILD_TESTING "Build the testing tree." ON)
include(CTest)
enable_testing()
message(STATUS "Testing is " ${BUILD_TESTING})
add_subdirectory(unittests)


#
# Add main source subdirectory
#
add_subdirectory(src)

add_subdirectory(include)

# Find MPI library
find_package(MPI REQUIRED)
target_link_libraries(lsmscore PUBLIC MPI::MPI_CXX)

# Find OpenMP library
option(BUILD_WITH_OPENMP "Build with OpenMP" ON)

if(BUILD_WITH_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(lsmscore PUBLIC OpenMP::OpenMP_CXX)
endif()

# Find HDF5
set(HDF5_PREFER_PARALLEL false)
find_package(HDF5 REQUIRED)
target_link_libraries(lsmscore PUBLIC HDF5::HDF5)

# Lapack
if (NOT DEFINED BLA_VENDOR)
    set(BLA_VENDOR "Generic")
endif ()
find_package(LAPACK REQUIRED)
target_link_libraries(lsmscore PUBLIC LAPACK::LAPACK)

# BLAS
find_package(BLAS REQUIRED)
target_link_libraries(lsmscore PUBLIC BLAS::BLAS)

# Lua
target_link_libraries(lsmscore PUBLIC Lua::Lua)

# mjson
target_link_libraries(lsmscore PUBLIC mjson)

# LibXC
# target_link_libraries(lsmscore PUBLIC libxc::libxc)

target_include_directories(lsmscore PUBLIC ${libxc_INCLUDE_DIR})
target_link_libraries(lsmscore PUBLIC ${libxc_LIBRARIES})
target_compile_definitions(lsmscore PUBLIC USE_LIBXC)

# DL
target_link_libraries(lsmscore PUBLIC ${CMAKE_DL_LIBS})
target_compile_definitions(lsmscore PUBLIC HAS_BACKTRACE)

# Include the common files
target_include_directories(lsmscore PUBLIC ${CMAKE_SOURCE_DIR}/include)


