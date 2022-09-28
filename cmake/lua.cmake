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

target_compile_definitions(lua
        PRIVATE
        $<$<PLATFORM_ID:Linux>:LUA_USE_LINUX LUA_COMPAT_5_2>)

target_compile_definitions(lua
        PRIVATE
        $<$<PLATFORM_ID:APPLE>:LUA_USE_MACOSX>)

target_compile_options(lua
        PRIVATE
        $<$<OR:$<C_COMPILER_ID:AppleClang>,$<C_COMPILER_ID:Clang>,$<C_COMPILER_ID:GNU>>:
        -Wextra -Wshadow -Wsign-compare -Wundef -Wwrite-strings -Wredundant-decls
        -Wdisabled-optimization -Waggregate-return -Wdouble-promotion -Wdeclaration-after-statement
        -Wmissing-prototypes -Wnested-externs -Wstrict-prototypes -Wc++-compat -Wold-style-definition>)

target_include_directories(lua
        INTERFACE
        ${srcDir})

add_library(Lua::Lua ALIAS lua)

