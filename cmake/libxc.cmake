#
# Created by Franco P. Moitzi
#

#
# Libxc libaries
#

if (DEFINED libxc_LIBRARIES)
    get_filename_component(libxc_PARENT_PATH ${libxc_LIBRARIES} DIRECTORY )
endif()


set(libxc_FOUND true)

find_library(libxc_LIBRARIES NAMES xc PATHS ${CMAKE_BINARY_DIR}/external/libxc/lib/ ${libxc_PARENT_PATH})

if (NOT libxc_LIBRARIES)
    message(STATUS "Libxc library was not found")
    set(libxc_FOUND false)
    set(libxc_LIBRARIES
            ${CMAKE_BINARY_DIR}/external/libxc/lib/${CMAKE_STATIC_LIBRARY_PREFIX}xc${CMAKE_STATIC_LIBRARY_SUFFIX}
            CACHE FILEPATH "libxc file path" FORCE)
endif()

find_path(libxc_INCLUDE_DIR NAMES xc.h PATHS ${CMAKE_BINARY_DIR}/external/libxc/include ${libxc_INCLUDE_DIR})

if (NOT libxc_INCLUDE_DIR)
    message(STATUS "Libxc include dir was not found")
    set(libxc_FOUND false)
    set(libxc_INCLUDE_DIR
            ${CMAKE_BINARY_DIR}/external/libxc/include CACHE PATH "libxc include dir" FORCE)
endif()

if (EXISTS ${libxc_LIBRARIES} AND EXISTS ${libxc_INCLUDE_DIR})
    set(libxc_FOUND true)
    message(STATUS "libxc was found")
else()
    set(libxc_FOUND false)
    message(STATUS "libxc was not found and will be installed")
endif ()

message(STATUS "libxc library: " ${libxc_LIBRARIES})
message(STATUS "libxc include: " ${libxc_INCLUDE_DIR})

if (NOT libxc_FOUND)

    cmake_policy(SET CMP0111 NEW)

    find_program(AUTORECONF_EXECUTABLE
            NAMES autoreconf
            DOC "Autoreconf" REQUIRED)

    find_program(AUTOCONF_EXECUTABLE
            NAMES autoconf
            DOC "Autoconf" REQUIRED)

    find_program(AUTOMAKE_EXECUTABLE
            NAMES automake
            DOC "Automake" REQUIRED)

    find_program(MAKE_EXECUTABLE
            NAMES gmake make
            NAMES_PER_DIR
            DOC "GNU Make")

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(Autotools
            REQUIRED_VARS AUTOCONF_EXECUTABLE AUTOMAKE_EXECUTABLE MAKE_EXECUTABLE)

    file(COPY ${PROJECT_SOURCE_DIR}/external/libxc-6.1.0
            DESTINATION ${CMAKE_BINARY_DIR}/external)

    set(_src ${CMAKE_BINARY_DIR}/external/libxc-6.1.0)
    get_filename_component(_src "${_src}" REALPATH)

    set(_install ${CMAKE_BINARY_DIR}/external/libxc)
    file(MAKE_DIRECTORY ${_install})
    file(MAKE_DIRECTORY ${_install}/include)
    get_filename_component(_install "${_install}" REALPATH)
    include(ExternalProject)

    ExternalProject_Add(libxc
            SOURCE_DIR ${_src}
            BUILD_IN_SOURCE true
            CONFIGURE_COMMAND ${AUTORECONF_EXECUTABLE} -i
            COMMAND ./configure --prefix=${_install} CC=${CMAKE_C_COMPILER}
            BUILD_COMMAND ${MAKE_EXECUTABLE}
            INSTALL_COMMAND ${MAKE_EXECUTABLE} install
            BUILD_BYPRODUCTS ${_install}/lib/libxc.a ${_install}/lib/libxc.so
            )


endif ()


add_library(libxc::libxc SHARED IMPORTED GLOBAL)
set_target_properties(libxc::libxc PROPERTIES IMPORTED_LOCATION ${libxc_LIBRARIES})
target_include_directories(libxc::libxc INTERFACE ${libxc_INCLUDE_DIR})
target_compile_definitions(libxc::libxc INTERFACE USE_LIBXC)
add_dependencies(libxc::libxc libxc)
