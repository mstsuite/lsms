
#ifndef LSMS_READ_INPUT_H
#define LSMS_READ_INPUT_H

#include "LSMSMode.hpp"
#include "LuaInterface/LuaSupport.hpp"
#include "MixerType.hpp"
#include "MixingParameter.hpp"
#include "Potential/PotentialShifter.hpp"
#include "SystemParameters.hpp"
#include "mixing_params.hpp"

int readInput(lua_State *L,
              LSMSSystemParameters &lsms,
              CrystalParameters &crystal,
              lsms::MixingParameterPack &mix,
              PotentialShifter &potentialShifter,
              AlloyMixingDesc &alloyDesc);

int readInput(lua_State *L,
              LSMSSystemParameters &lsms,
              CrystalParameters &crystal,
              MixingParameters &mix,
              PotentialShifter &potentialShifter,
              AlloyMixingDesc &alloyDesc);

#endif
