

find_program(MAKE_EXECUTABLE NAMES gmake make REQUIRED)

include(ExternalProject)


set(Lua_LIBRARIES
        ${PROJECT_BINARY_DIR}/external/lua/lib/${CMAKE_STATIC_LIBRARY_PREFIX}lua${CMAKE_STATIC_LIBRARY_SUFFIX})

set(Lua_INCLUDE_DIR
        ${PROJECT_BINARY_DIR}/external/lua/include
        )


if (EXISTS ${Lua_LIBRARIES} AND EXISTS ${Lua_INCLUDE_DIR})
    set(Lua_FOUND true)
    message(STATUS "Lua library was already installed")
endif ()

if (NOT Lua_FOUND)

    file(COPY ${PROJECT_SOURCE_DIR}/external/lua-5.2.4
            DESTINATION ${PROJECT_BINARY_DIR}/external)

    set(_src ${PROJECT_BINARY_DIR}/external/lua-5.2.4)
    get_filename_component(_src "${_src}" REALPATH)

    set(_install ${PROJECT_BINARY_DIR}/external/lua)
    file(MAKE_DIRECTORY ${_install})
    get_filename_component(_install "${_install}" REALPATH)

    ExternalProject_Add(Lua
            SOURCE_DIR ${_src}
            BUILD_IN_SOURCE false
            CONFIGURE_COMMAND sed -i "/^INSTALL_TOP/c INSTALL_TOP=${_install}" ${_src}/Makefile
            BUILD_COMMAND ${MAKE_EXECUTABLE} -C ${_src}
            INSTALL_COMMAND ${MAKE_EXECUTABLE} install -C ${_src}
            )



endif ()

add_library(Lua::Lua INTERFACE IMPORTED GLOBAL)
target_include_directories(Lua::Lua INTERFACE ${Lua_INCLUDE_DIR})
target_link_libraries(Lua::Lua INTERFACE ${Lua_LIBRARIES})

if (NOT Lua_FOUND)
    add_dependencies(Lua::Lua Lua)
endif ()