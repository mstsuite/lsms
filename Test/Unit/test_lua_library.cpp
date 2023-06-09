//
// Created by F.Moitzi on 18.09.2021.
//

#include <gtest/gtest.h>

#include <cstring>
#include <filesystem>
#include <iostream>
#include <lua.hpp>

#include "MixerType.hpp"
#include "MixingParameter.hpp"
#include "accel_common.hpp"

namespace lua_libarary_tests {

TEST(LuaLibrary, Test1) {
  // initialization
  lua_State *L = luaL_newstate();
  luaL_openlibs(L);

  // execute script
  const char lua_script[] =
      "function sum(a, b) return a+b; end";  // a function that returns sum of
  // two
  int load_stat =
      luaL_loadbuffer(L, lua_script, strlen(lua_script), lua_script);
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

    if (sumval != 11.0) {
      FAIL() << "Sum isn't correct";
    }
  }

  // cleanup
  lua_close(L);
}

TEST(LuaLibrary, Test2) {
  std::string inputFileName("input.lua");

  ASSERT_TRUE(std::filesystem::exists(inputFileName));

  lua_State *L = luaL_newstate();
  luaL_openlibs(L);

  luaL_loadfile(L, inputFileName.c_str());
  lua_pcall(L, 0, 0, 0);

  lsms::MixingParameterPack mixingParams;

  lua_getglobal(L, "mixing");
  lua_istable(L, -1);
  std::size_t len = lua_rawlen(L, -1);

  for (auto index = 0; index < len; index++) {
    // Our actual index will be +1 because Lua 1 indexes tables.
    int actualIndex = index + 1;
    lua_pushinteger(L, actualIndex);
    lua_gettable(L, -2);

    // Quantity checks
    lua_pushstring(L, "quantity");
    lua_gettable(L, -2);

    std::string type = lua_tostring(L, -1);

    lua_pop(L, 1);

    // Mixing algorithm
    lua_pushstring(L, "algorithm");
    lua_gettable(L, -2);

    std::string algorithm = lua_tostring(L, -1);

    auto algo = lsms::getMixerType(algorithm);

    if (type == "potential") {
      mixingParams.pot_mixer_type = algo;
    } else if (type == "density") {
      mixingParams.chd_mixer_type = algo;
    } else {
      throw std::runtime_error("No implemented");
    }

    lua_pop(L, 1);

    // Mixing parameters
    lua_pushstring(L, "alpha");
    lua_gettable(L, -2);
    if (lua_type(L, -1) != LUA_TNIL) {
      double alpha = lua_tonumber(L, -1);

      if (type == "potential") {
        mixingParams.potMixingParameter.alpha = alpha;
      } else if (type == "density") {
        mixingParams.chdMixingParameter.alpha = alpha;
      } else {
        throw std::runtime_error("No implemented");
      }
    };
    lua_pop(L, 1);

    lua_pushstring(L, "w0");
    lua_gettable(L, -2);
    if (lua_type(L, -1) != LUA_TNIL) {
      double w0 = lua_tonumber(L, -1);

      if (type == "potential") {
        mixingParams.potMixingParameter.w0 = w0;
      } else if (type == "density") {
        mixingParams.chdMixingParameter.w0 = w0;
      } else {
        throw std::runtime_error("No implemented");
      }
    };
    lua_pop(L, 1);

    lua_pushstring(L, "max_broyden");
    lua_gettable(L, -2);
    if (lua_type(L, -1) != LUA_TNIL) {
      unsigned int max_broyden = lua_tonumber(L, -1);

      if (type == "potential") {
        mixingParams.potMixingParameter.max_broyden = max_broyden;
      } else if (type == "density") {
        mixingParams.chdMixingParameter.max_broyden = max_broyden;
      } else {
        throw std::runtime_error("No implemented");
      }
    }
    lua_pop(L, 1);

    lua_pushstring(L, "iter_reset");
    lua_gettable(L, -2);
    if (lua_type(L, -1) != LUA_TNIL) {
      unsigned int iter_reset = lua_tonumber(L, -1);

      if (type == "potential") {
        mixingParams.potMixingParameter.iter_reset = iter_reset;
      } else if (type == "density") {
        mixingParams.chdMixingParameter.iter_reset = iter_reset;
      } else {
        throw std::runtime_error("No implemented");
      }
    }
    lua_pop(L, 1);

    lua_pop(L, 1);
  }

  ASSERT_EQ(mixingParams.chd_mixer_type, lsms::MixerType::BROYDEN_MIXER);
  ASSERT_EQ(mixingParams.pot_mixer_type, lsms::MixerType::NO_MIXER);
  ASSERT_EQ(mixingParams.potMixingParameter.alpha, 0.54);
  ASSERT_EQ(mixingParams.chdMixingParameter.alpha, 0.03);
  ASSERT_EQ(mixingParams.chdMixingParameter.max_broyden, 8);

  lua_close(L);
}

}  // namespace lua_libarary_tests