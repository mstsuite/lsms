
#ifndef LSMS_READ_INPUT_H
#define LSMS_READ_INPUT_H

#include "SystemParameters.hpp"
#include "mixing.hpp"
#include "LSMSMode.hpp"
#include "LuaInterface/LuaSupport.hpp"
#include "Potential/PotentialShifter.hpp"
#include "mixing_params.hpp"

int readInput(lua_State *L,
              LSMSSystemParameters &lsms,
              CrystalParameters &crystal,
              MixingParameters &mix,
              PotentialShifter &potentialShifter,
              AlloyMixingDesc &alloyDesc);


#endif
