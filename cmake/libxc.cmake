#
# Created by Franco P. Moitzi
# Edited by Zhantong Qiu
#

# Check if libxc is already available
if (DEFINED libxc_LIBRARIES)
    get_filename_component(libxc_PARENT_PATH ${libxc_LIBRARIES} DIRECTORY)
endif()

set(libxc_BUILD_SRC_PATH ${CMAKE_BINARY_DIR}/external/libxc-6.1.0)
set(libxc_LIBRARIES_PATH ${CMAKE_BINARY_DIR}/external/libxc-build)
set(libxc_INCLUDE_DIR_PATH ${CMAKE_BINARY_DIR}/external/libxc-build/src)

# check if neccessary files are available
find_library(libxc_LIBRARIES NAMES xc PATHS ${libxc_LIBRARIES_PATH} ${libxc_PARENT_PATH})
find_path(libxc_INCLUDE_DIR NAMES xc.h PATHS ${libxc_INCLUDE_DIR_PATH} ${libxc_INCLUDE_DIR})
find_path(libxc_INCLUDE_BUILD_DIR NAMES xc_version.h PATHS ${libxc_LIBRARIES_PATH} ${libxc_INCLUDE_BUILD_DIR})

message(STATUS "libxc_LIBRARIES: ${libxc_LIBRARIES}")
message(STATUS "libxc_INCLUDE_DIR: ${libxc_INCLUDE_DIR}")
message(STATUS "libxc_INCLUDE_BUILD_DIR: ${libxc_INCLUDE_BUILD_DIR}")

if (EXISTS ${libxc_LIBRARIES} AND EXISTS ${libxc_INCLUDE_DIR} AND EXISTS ${libxc_INCLUDE_BUILD_DIR})
    message(STATUS "libxc was found")
    
    # Create imported target for existing installation
    add_library(libxc::libxc SHARED IMPORTED GLOBAL)
    set_target_properties(libxc::libxc PROPERTIES IMPORTED_LOCATION ${libxc_LIBRARIES})
    target_include_directories(libxc::libxc INTERFACE ${libxc_BUILD_SRC_PATH}/src ${libxc_INCLUDE_DIR} ${libxc_INCLUDE_BUILD_DIR})
        
    target_compile_definitions(libxc::libxc INTERFACE USE_LIBXC)
else()
    message(STATUS "libxc was not found and will be built from source")
    
    # Build libxc from source
    # Copy the source file to the build directory
    file(COPY ${PROJECT_SOURCE_DIR}/external/libxc-6.1.0
    DESTINATION ${CMAKE_BINARY_DIR}/external)
    
    set(_src ${libxc_BUILD_SRC_PATH})
    get_filename_component(_src "${_src}" REALPATH)
    
    # Add subdirectory to build the native target
    # This will include the CMakeLists.txt in the libxc source directory
    # So we can use the target libxc::xc
    add_subdirectory(${_src} ${libxc_LIBRARIES_PATH})
    
    # Set proper include directories for both source and generated files
    set(libxc_INCLUDE_DIRS
        ${_src}/src                 # Original source headers
        ${libxc_INCLUDE_DIR_PATH}   # Generated xc.h
        ${libxc_LIBRARIES_PATH}     # Generated xc_version.h
    )
    
    # Create alias target and set include directories
    add_library(libxc::libxc ALIAS xc)
    
    target_include_directories(xc PUBLIC 
        $<BUILD_INTERFACE:${_src}/src>
        $<BUILD_INTERFACE:${libxc_INCLUDE_DIR_PATH}>
        $<BUILD_INTERFACE:${libxc_LIBRARIES_PATH}>)
endif()

