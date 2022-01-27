//
// Created by F.Moitzi on 18.09.2021.
//

#include <gtest/gtest.h>

#include <cstring>
#include <cstdlib>

#include <array>
#include <vector>

#include <lua.hpp>

TEST(LuaTest, LuaScripting1) {

  // initialization
  lua_State *L = luaL_newstate();
  luaL_openlibs(L);

  // execute script
  const char lua_script[] = "function sum(a, b) return a+b; end"; // a function that returns sum of two
  int load_stat = luaL_loadbuffer(L, lua_script, strlen(lua_script), lua_script);
  lua_pcall(L, 0, 0, 0);

  // load the function from global
  lua_getglobal(L, "sum");
  if (lua_isfunction(L, -1)) {
    // push function arguments into stack
    lua_pushnumber(L, 5.0);
    lua_pushnumber(L, 6.0);
    lua_pcall(L, 2, 1, 0);
    double sumval = 0.0;
    if (!lua_isnil(L, -1)) {
      sumval = lua_tonumber(L, -1);
      lua_pop(L, 1);
    }

    ASSERT_EQ(sumval, 11.0);


  }

  lua_close(L);
}