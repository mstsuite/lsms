# Source Lists -----------------------------------------------------------------

set(LSMS_LUA_INTERFACE_SOURCES
    LuaInterface.cpp
    LuaSupport.cpp
)

# Libraries --------------------------------------------------------------------
add_library(LuaInterface STATIC ${LSMS_LUA_INTERFACE_SOURCES})

target_link_libraries(LuaInterface PUBLIC Lua)

target_include_directories(LuaInterface PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)
