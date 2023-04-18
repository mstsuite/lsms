
#ifndef LSMS_RADIALGRID_LUA_HPP
#define LSMS_RADIALGRID_LUA_HPP

#include "RadialGrid/RadialGrid.hpp"
#include "lauxlib.h"
#include "lua.h"
#include "lualib.h"

// a radial grid in lua is represented as a table that contains the radial grid
// pointer as light userdata and the appropriate metatable "LSMS.RadialGrid".
// radial grids that have been allocated by lua are considered temporary and
// contain the field "ownedByLua" set to true. only LSMS internal functions are
// allowed to assign pointer values and care should be taken that pointers owned
// by lua are never passed to such functions!

struct {
  RadialGrid *g;
  bool ownedByLua;
} RadialGridLua;

static const struct luaL_Reg RadialGrid_lib_f[] = {
    {"new", newTemporaryRadialGridLua}, {NULL, NULL}};

static const struct luaL_Reg RadialGrid_lib_m[] = {
    {"__index", radialGridLua__index},
    {"generate", generateRadialGridLua},
    {"__tostring", luaRadialGridToString},
    {"__gc", radialGridLua__gc},
    {"size", sizeRadialGridLua},
    {"get", getRadialGridLua},
    {"x", getRadialGridXLua},
    //  {"copy", copyRadialGridLua},
    {NULL, NULL}};

int newTemporaryRadialGridLua(lua_State *L);

int radialGridLua__gc(lua_State *L)

    int generateRadialGridLua(lua_State *L);

int sizeRadialGridLua(lua_State *L);

int getRadialGridLua(lua_State *L);

int getRadialGridXLua(lua_State *L);

int luaRadialGridToString(lua_State *L);

int radialGridLua__index(lua_State *L);

int luaopen_RadialGrid(lua_State *L);

#endif  // LSMS_RADIALGRID_LUA_HPP