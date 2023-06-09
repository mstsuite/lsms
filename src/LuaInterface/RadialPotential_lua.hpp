
#ifndef LSMS_RADIALPOTENTIAL_LUA_HPP
#define LSMS_RADIALPOTENTIAL_LUA_HPP

#include "RadialGrid/RadialGrid.hpp"
#include "RadialGrid/RadialPotential.hpp"
#include "lauxlib.h"
#include "lua.h"
#include "lualib.h"

// a radial grid in lua is represented as a table that contains the radial grid
// pointer as light userdata and the appropriate metatable "LSMS.RadialGrid".
// radial grids that have been allocated by lua are considered temporary and
// contain the field "ownedByLua" set to true. only LSMS internal functions are
// allowed to assign pointer values and care should be taken that pointers owned
// by lua are never passed to such functions!

int newTemporaryRadialPotentialLua(lua_State *L);

int radialPotentialLua__gc(lua_State *L);

int getRadialGridFromPotentialLua(lua_State *L);

int getRadialPotentialLua(lua_State *L);

int syncRadialPotentialLua(lua_State *L);

int sizeRadialPotentialLua(lua_State *L);

int luaRadialPotentialToString(lua_State *L);

int radialPotentialLua__index(lua_State *L);

int luaopen_RadialPotential(lua_State *L);

static const struct luaL_Reg RadialPotential_lib_f[] = {
    {"new", newTemporaryRadialPotentialLua}, {NULL, NULL}};

static const struct luaL_Reg RadialPotential_lib_m[] = {
    {"__index", radialPotentialLua__index},
    //  {"generate", generateRadialGridLua},
    {"__tostring", luaRadialPotentialToString},
    {"__gc", radialPotentialLua__gc},
    {"grid", getRadialGridFromPotentialLua},
    {"size", sizeRadialPotentialLua},
    {"get", getRadialPotentialLua},
    {"sync", syncRadialPotentialLua},
    //  {"copy", copyRadialGridLua},
    {NULL, NULL}};

#endif  // LSMS_RADIALPOTENTIAL_LUA_HPP