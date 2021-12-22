

if (DEFINED Lua_LIBRARIES AND DEFINED Lua_INCLUDE_DIR)
    if (EXISTS ${Lua_LIBRARIES} AND EXISTS ${Lua_INCLUDE_DIR})
        message(STATUS "Lua path was correctly defined")
    else()
        message(ERROR "Specified path for lua libaries is wrong")
    endif ()
endif()

if (NOT DEFINED Lua_LIBRARIES)
    set(Lua_LIBRARIES
            ${CMAKE_BINARY_DIR}/external/lua/lib/${CMAKE_STATIC_LIBRARY_PREFIX}lua${CMAKE_STATIC_LIBRARY_SUFFIX})

    set(Lua_INCLUDE_DIR
            ${CMAKE_BINARY_DIR}/external/lua/include
            )
endif ()

if (EXISTS ${Lua_LIBRARIES} AND EXISTS ${Lua_INCLUDE_DIR})
    set(Lua_FOUND true)
    message(STATUS "Lua was found")
    message(STATUS "Lua library: " ${Lua_LIBRARIES})
    message(STATUS "Lua include: " ${Lua_INCLUDE_DIR})
else()
    set(Lua_FOUND false)
    message(STATUS "Lua was not found and will be installed")
endif ()

if (NOT Lua_FOUND)
    find_program(MAKE_EXECUTABLE NAMES gmake make REQUIRED)

    file(COPY ${PROJECT_SOURCE_DIR}/external/lua-5.2.4
            DESTINATION ${CMAKE_BINARY_DIR}/external)

    set(_src ${CMAKE_BINARY_DIR}/external/lua-5.2.4)
    get_filename_component(_src "${_src}" REALPATH)

    set(_install ${CMAKE_BINARY_DIR}/external/lua)
    file(MAKE_DIRECTORY ${_install})
    get_filename_component(_install "${_install}" REALPATH)

    set(_include ${CMAKE_BINARY_DIR}/external/lua/include)
    file(MAKE_DIRECTORY ${_include})
    get_filename_component(_include "${_include}" REALPATH)

    include(ExternalProject)

    ExternalProject_Add(Lua
            SOURCE_DIR ${_src}
            BUILD_IN_SOURCE true
            CONFIGURE_COMMAND sed -i '' -e "/^INSTALL_TOP/c INSTALL_TOP=${_install}" ${_src}/Makefile
            BUILD_COMMAND ${MAKE_EXECUTABLE} -C ${_src}
            INSTALL_COMMAND ${MAKE_EXECUTABLE} install -C ${_src}
            )
endif()


# convert files outside of a CMake project into logical targets inside of the project.
# No build files are created. Global makes it visible every where
add_library(Lua::Lua STATIC IMPORTED GLOBAL)
set_target_properties(Lua::Lua PROPERTIES IMPORTED_LOCATION ${Lua_LIBRARIES})
target_include_directories(Lua::Lua INTERFACE ${Lua_INCLUDE_DIR})
add_dependencies(Lua::Lua Lua)