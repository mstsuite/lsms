# prohibit in-source-builds
if (${CMAKE_BINARY_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
	message(STATUS "In-source-builds are not allowed")
	message(STATUS "Clean your source directory (e.g. delete the CMakeCache.txt file)")
	message(FATAL_ERROR "Please create a separate build directory and call CMake again")
endif()

cmake_minimum_required(VERSION 3.12)



if (WIN32 AND (CMAKE_GENERATOR STREQUAL "Ninja"))
	set(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE})
endif()

# Variables --------------------------------------------------------------------

set(PROJECT_NAME LSMS)

# Project ----------------------------------------------------------------------
project(${PROJECT_NAME})

enable_language (Fortran)

# Set a default build type if none was specified
set(default_build_type "Release")
if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
  set(default_build_type "Debug")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

#RPATH
if(UNIX AND NOT APPLE)
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    set(CMAKE_INSTALL_RPATH "$ORIGIN/lib/:$$ORIGIN/lib/")
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 99)

# Options ----------------------------------------------------------------------
option(ENABLE_OPENMP "Enable OpenMP" OFF)
option(ENABLE_CUDA "Enable CUDA" OFF)

# Dependencies -----------------------------------------------------------------
find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

if(ENABLE_OPENMP)
    find_package(OpenMP)
endif()

if(ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
    else()
        message(FATAL_ERROR "Could not find CUDA compiler")
    endif()
endif()

# Lua
find_package(Lua QUIET)
if(NOT LUA_FOUND)
    message(STATUS "Could not find Lua library")
    message(STATUS "\tUsing build in Lua 5.3.5")
    include(ExternalProject)
    find_program(MAKE_EXE NAMES gmake nmake make)
    ExternalProject_Add(LUA
        URL https://www.lua.org/ftp/lua-5.3.5.tar.gz
        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/lua-5.3.5
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ${MAKE_EXE} posix
        BUILD_IN_SOURCE True
        INSTALL_COMMAND ${MAKE_EXE} install INSTALL_TOP=${CMAKE_CURRENT_BINARY_DIR}/lua-5.3.5
    )

    set(LUA_LIBRARIES liblua.a)
    set(LUA_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/lua-5.3.5/include)

    # make sure the include directory exists
    file(MAKE_DIRECTORY ${LUA_INCLUDE_DIR})

    add_library(Lua STATIC IMPORTED)
    set_target_properties(Lua PROPERTIES
        IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/lua-5.3.5/lib/liblua.a
        INTERFACE_INCLUDE_DIRECTORIES ${LUA_INCLUDE_DIR}
    )
endif()

# Settings ---------------------------------------------------------------------
set(CMAKE_VERBOSE_MAKEFILE ON)

# Subdirectories ---------------------------------------------------------------
add_subdirectory(mjson)
add_subdirectory(src/RadialGrid)
add_subdirectory(src/VORPOL)
add_subdirectory(src/Madelung)
add_subdirectory(src/Potential)
add_subdirectory(src/SingleSite)
add_subdirectory(src/MultipleScattering)
add_subdirectory(src/Core)
add_subdirectory(src/LuaInterface)
add_subdirectory(src/Communication)
add_subdirectory(src/Misc)
add_subdirectory(src/Main)
add_subdirectory(src/Accelerator)
add_subdirectory(src/TotalEnergy)


# Libraries --------------------------------------------------------------------

# Executables ------------------------------------------------------------------

