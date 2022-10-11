//
// Created by F.Moitzi on 21.12.2022.
//

#ifndef LSMS_MIXINGPARAMETER_HPP
#define LSMS_MIXINGPARAMETER_HPP

#include "MixerType.hpp"
#include "utils.h"

namespace lsms {

constexpr static unsigned int DEFAULT_INIT_MIXER_TYPE = MixerType::SIMPLE_MIXER;

constexpr static unsigned int DEFAULT_CHD_MIXER_TYPE = MixerType::SIMPLE_MIXER;

constexpr static unsigned int DEFAULT_SPD_MIXER_TYPE = MixerType::SIMPLE_MIXER;

constexpr static unsigned int DEFAULT_INIT_SPD_MIXER_TYPE = MixerType::SIMPLE_MIXER;

constexpr static unsigned int DEFAULT_POT_MIXER_TYPE = MixerType::BROYDEN_MIXER;

class MixingParameter {
 public:
  double alpha;
  double w0;
  unsigned int max_broyden;
  unsigned int iter_reset;

  MixingParameter() = default;

  explicit MixingParameter(double alpha,
                           double w0 = 0.0,
                           unsigned int max_broyden = 0,
                           unsigned int iter_reset = 0
  );
};

struct MixingParameterPack {

  // Number of initial spin iterations
  unsigned int n_init_spin_iterations = 0;

  // Number of initial iterations
  unsigned int n_init_iterations = 10;

  // Type of initial mixer (spd and chd)
  unsigned int init_mixer_type = DEFAULT_INIT_MIXER_TYPE;

  // Initial mixer parameters
  MixingParameter initMixingParameter{1.0, 0.01, 10, 25};

  // Type of charge density mixer
  unsigned int chd_mixer_type = DEFAULT_CHD_MIXER_TYPE;

  // Charge density mixing parameters
  MixingParameter chdMixingParameter{1.0, 0.01, 10, 25};

  // Type of spin density mixer
  unsigned int spd_mixer_type = DEFAULT_SPD_MIXER_TYPE;

  // Spin density mixing parameters
  MixingParameter spdMixingParameter{1.0, 0.01, 10, 25};

  // Type of spin density mixer
  unsigned int init_spd_mixer_type = DEFAULT_INIT_SPD_MIXER_TYPE;

  // Spin density mixing parameters
  MixingParameter initSpdMixingParameter{1.0, 0.01, 10, 25};

  // Type of potential mixer
  unsigned int pot_mixer_type = DEFAULT_POT_MIXER_TYPE;

  // Potential mixing parameters
  MixingParameter potMixingParameter{0.05, 0.01, 10, 25};
};

}  // namespace lsms

#endif  // LSMS_MIXINGPARAMETER_HPP
