//
// Created by F.Moitzi on 21.12.2022.
//

#ifndef LSMS_MIXERTYPE_H
#define LSMS_MIXERTYPE_H

#include <string>

namespace lsms {

struct MixerType {
  static constexpr const unsigned int NO_MIXER = 0;
  static constexpr const unsigned int SIMPLE_MIXER = 1;
  static constexpr const unsigned int BROYDEN_MIXER = 2;
};

std::string getMixerName(unsigned int mixerType);

unsigned int getMixerType(std::string &name);

}  // namespace lsms

#endif  // LSMS_MIXERTYPE_H
