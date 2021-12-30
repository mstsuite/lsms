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

# Enable compilation for Accelerators
option(USE_ACCELERATOR_CUDA_C "Use the cuda accelerator" OFF)
# option(USE_BUILDKKRMATRIX_GPU "Build also the kkr matrix on GPU" ON)
option(USE_ACCELERATOR_HIP "Use HIP for accelerator" OFF)

if(USE_ACCELERATOR_CUDA_C AND USE_ACCELERATOR_HIP)
  message(FATAL_ERROR "Both ACCELERATOR_CUDA_C and ACCELERATOR_HIP defined.\nOnly ONE of these allowed!")
endif()

if (USE_ACCELERATOR_CUDA_C)
    message(STATUS "ACCELERATOR_CUDA_C was activated")
    enable_language(CUDA)
endif ()

option(USE_ACCELERATOR_HIP "Use HIP for accelerator" OFF)

if (USE_ACCELERATOR_HIP)
    message(STATUS "ACCELERATOR_HIP was activated")
endif ()



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
option(USE_LIBXC "Use the libxc libary" ON)
if (USE_LIBXC)
    message(STATUS "LIBXC is used")
    include(libxc)
else()
    message(STATUS "LIBXC is not used")
endif ()

# Build mjson library
include(mjson)


#
# Enable Testing with CTest
#
option(BUILD_TESTING "Build the testing tree." OFF)
if (BUILD_TESTING)
    include(CTest)
    enable_testing()
    message(STATUS "Testing is " ${BUILD_TESTING})
    add_subdirectory(unittests)
endif ()

#
# Add main source subdirectory
#
add_subdirectory(src)

add_subdirectory(include)

#
#
#


# Find OpenMP library
option(BUILD_WITH_OPENMP "Build with OpenMP" ON)
if (BUILD_WITH_OPENMP)
    find_package(OpenMP REQUIRED)
    target_link_libraries(lsmscore PUBLIC OpenMP::OpenMP_CXX)
endif ()

#
# CUDA
#
if (USE_ACCELERATOR_CUDA_C)

    target_compile_definitions(lsmscore PUBLIC ACCELERATOR_CUDA_C)

    option(CUDA_DEBUG "Debug CUDA routines" OFF)

    if (CUDA_DEBUG)
        target_compile_definitions(lsmscore PUBLIC "CUDA_DEBUG")
    endif ()

    #    if (USE_BUILDKKRMATRIX_GPU)
    #        message(STATUS "KKR Matrix build on GPU")
    #        target_compile_definitions(lsmscore PUBLIC BUILDKKRMATRIX_GPU)
    #    else ()
    #        message(STATUS "KKR Matrix is not build on GPU")
    #    endif ()

    if (BUILD_WITH_OPENMP)
        message(STATUS "OpenMP is also defined for GPU")
        target_compile_options(lsmscore PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler -fopenmp>)
    endif ()

    if (NOT DEFINED CMAKE_CUDA_STANDARD)
        set(CMAKE_CUDA_STANDARD 14)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    endif ()

    if (NOT DEFINED CMAKE_CUDA_SEPARABLE_COMPILATION)
        set(CMAKE_CUDA_SEPARABLE_COMPILATION OFF)
    endif ()

    find_package(CUDAToolkit REQUIRED)

    message(STATUS "CUDA Version: " ${CUDAToolkit_VERSION})
    message(STATUS "CUDA Include dirs: " ${CUDAToolkit_INCLUDE_DIRS})
    message(STATUS "CUDA Library: " ${CUDAToolkit_LIBRARY_DIR})

    target_link_libraries(lsmscore PUBLIC CUDA::toolkit)

    # Check if the necessary libraries have been found
    if (TARGET CUDA::cublas)
        target_link_libraries(lsmscore PUBLIC CUDA::cublas)
    else ()
        message(FATAL_ERROR "CuBLAS was not defined or found")
    endif ()

    if (TARGET CUDA::cusolver)
        target_link_libraries(lsmscore PUBLIC CUDA::cusolver)
    else ()
        message(FATAL_ERROR "Cusolver was not defined or found")
    endif ()
endif ()

if (USE_ACCELERATOR_HIP)

    target_compile_definitions(lsmscore PUBLIC ACCELERATOR_HIP)

endif ()

if (DEFINED MST_LINEAR_SOLVER_DEFAULT)
    message(STATUS "Default linear solver: " ${MST_LINEAR_SOLVER_DEFAULT})
    set(MST_LINEAR_SOLVER_DEFAULT ${MST_LINEAR_SOLVER_DEFAULT} CACHE STRING "Default linear solver")
    target_compile_definitions(lsmscore PUBLIC "MST_LINEAR_SOLVER_DEFAULT=${MST_LINEAR_SOLVER_DEFAULT}")
endif ()

if (DEFINED MST_BUILD_KKR_MATRIX_DEFAULT)
    message(STATUS "Default kkr matrix build: " ${MST_BUILD_KKR_MATRIX_DEFAULT})
    set(MST_BUILD_KKR_MATRIX_DEFAULT ${MST_BUILD_KKR_MATRIX_DEFAULT} CACHE STRING "Default linear solver")
    target_compile_definitions(lsmscore PUBLIC "MST_BUILD_KKR_MATRIX_DEFAULT=${MST_BUILD_KKR_MATRIX_DEFAULT}")
endif ()


option(ARCH_IBM "Is IBM architecture" OFF)

if (ARCH_IBM)
   message(STATUS "Use the IBM architecture")
   target_compile_definitions(lsmscore PUBLIC "ARCH_IBM")
endif()

# Find MPI library
find_package(MPI REQUIRED COMPONENTS CXX OPTIONAL_COMPONENTS)
target_link_libraries(lsmscore PUBLIC MPI::MPI_CXX)


# Find HDF5
set(HDF5_PREFER_PARALLEL false)
find_package(HDF5 REQUIRED)
target_link_libraries(lsmscore PUBLIC HDF5::HDF5)

# Linear algebra libraries
option(USE_ESSL "Use IBM ESSL library" OFF)

if(USE_ESSL)
    find_package(ESSL REQUIRED)
    set(ESSL_WORKAROUND True)

    target_link_libraries(lsmscore PUBLIC ESSL)
else()
# Lapack
  option(SEARCH_LAPACK "Search for LAPACK" ON)
  if (NOT SEARCH_LAPACK)
      add_library(LAPACK::LAPACK STATIC IMPORTED GLOBAL)
      set_target_properties(LAPACK::LAPACK PROPERTIES IMPORTED_LOCATION ${LAPACK_LIBRARIES})
      if (DEFINED LAPACK_LINKER_FLAGS)
          target_compile_options(LAPACK::LAPACK INTERFACE ${LAPACK_LINKER_FLAGS})
      endif ()
      if (DEFINED LAPACK_INCLUDE_DIR)
          target_include_directories(LAPACK::LAPACK INTERFACE ${LAPACK_INCLUDE_DIR})
      endif ()
else ()
      if (NOT DEFINED BLA_VENDOR)
          set(BLA_VENDOR "Generic")
      endif ()
      find_package(LAPACK REQUIRED)
  endif ()
  target_link_libraries(lsmscore PUBLIC LAPACK::LAPACK)

# BLAS
  option(SEARCH_BLAS "Search for BLAS" ON)
  if (NOT SEARCH_BLAS)
      add_library(BLAS::BLAS STATIC IMPORTED GLOBAL)
      set_target_properties(BLAS::BLAS PROPERTIES IMPORTED_LOCATION ${BLAS_LIBRARIES})
      if (DEFINED BLAS_LINKER_FLAGS)
          target_compile_options(BLAS::BLAS INTERFACE ${BLAS_LINKER_FLAGS})
      endif ()
      if (DEFINED BLAS_INCLUDE_DIR)
          target_include_directories(BLAS::BLAS INTERFACE ${BLAS_INCLUDE_DIR})
      endif ()
  else ()
      if (NOT DEFINED BLA_VENDOR)
          set(BLA_VENDOR "Generic")
      endif ()
      find_package(BLAS REQUIRED)
  endif ()
  target_link_libraries(lsmscore PUBLIC BLAS::BLAS)
endif()

# Lua
target_link_libraries(lsmscore PUBLIC Lua::Lua)

# mjson
target_link_libraries(lsmscore PUBLIC mjson)

# LibXC
if (USE_LIBXC)
    target_link_libraries(lsmscore PUBLIC libxc::libxc)
    target_compile_definitions(lsmscore PUBLIC USE_LIBXC)
endif ()

# DL
target_link_libraries(lsmscore PUBLIC ${CMAKE_DL_LIBS})
target_compile_definitions(lsmscore PUBLIC HAS_BACKTRACE)

# Include the common files
target_include_directories(lsmscore PUBLIC ${CMAKE_SOURCE_DIR}/include)


