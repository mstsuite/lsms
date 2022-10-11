//
// Created by F.Moitzi on 23.12.2022.
//

#include "MixerType.hpp"

#include <exception>
#include <string>

namespace lsms {

unsigned int getMixerType(std::string &name) {
  if (name == "broyden") {
    return MixerType::BROYDEN_MIXER;
  } else if (name == "simple") {
    return MixerType::SIMPLE_MIXER;
  } else {
    return MixerType::NO_MIXER;
  }
}

std::string getMixerName(unsigned int mixerType) {
  if (mixerType == MixerType::SIMPLE_MIXER) {
    return {"Simple linear mixer"};
  } else if (mixerType == MixerType::BROYDEN_MIXER) {
    return {"Broyden mixer"};
  } else if (mixerType == MixerType::NO_MIXER) {
    return {"No mixer"};
  } else {
    return {""};
  }
};

}  // namespace lsms
