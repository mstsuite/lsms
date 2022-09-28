#
# Created by Franco P. Moitzi
#

set(srcDir ${PROJECT_SOURCE_DIR}/external/lua-5.2.4/src)

set(srcFiles lapi.c lcode.c lctype.c ldebug.c ldo.c ldump.c lfunc.c lgc.c llex.c
        lmem.c lobject.c lopcodes.c lparser.c lstate.c lstring.c ltable.c
        ltm.c
        lundump.c
        lvm.c
        lzio.c
        lauxlib.c
        lbaselib.c
        lbitlib.c
        lcorolib.c
        ldblib.c
        liolib.c
        lmathlib.c
        loslib.c
        lstrlib.c
        ltablib.c
        loadlib.c
        linit.c
        lapi.h
        lauxlib.h
        lcode.h
        lctype.h
        ldebug.h
        ldo.h
        lfunc.h
        lgc.h
        llex.h
        llimits.h
        lmem.h
        lobject.h
        lopcodes.h
        lparser.h
        lstate.h
        lstring.h
        ltable.h
        ltm.h
        luaconf.h
        lua.h
        lualib.h
        lundump.h
        lvm.h
        lzio.h)

list(TRANSFORM srcFiles PREPEND ${srcDir}/)

add_library(lua STATIC)

target_sources(lua PUBLIC
        ${srcFiles}
        )

target_include_directories(lua
        PUBLIC    
        ${srcDir})


if (CMAKE_SYSTEM_NAME STREQUAL "Linux")

    message(STATUS "Lua: Use linux and POSIX")

    target_compile_definitions(lua PUBLIC LUA_USE_LINUX)
    target_link_libraries(lua PUBLIC ${CMAKE_DL_LIBS})

    find_package(Readline REQUIRED)
    message(STATUS "Readline library: " ${Readline_LIBRARY})
    message(STATUS "Readline include dir: " ${Readline_INCLUDE_DIR})

    target_link_libraries(lua PUBLIC ${Readline_LIBRARY})
    target_include_directories(lua PUBLIC ${Readline_INCLUDE_DIR})

endif ()


add_library(Lua::Lua ALIAS lua)

