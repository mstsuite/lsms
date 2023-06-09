//
// Created by F.Moitzi on 13.11.2022.
//

#ifndef LSMS_STATES_HPP
#define LSMS_STATES_HPP

#include <vector>

namespace lsms {

class States {
 public:
  static void relativistic_atomic_states(int core_charge, std::vector<int> &n,
                                         std::vector<int> &l,
                                         std::vector<int> &spin,
                                         std::vector<int> &kappa,
                                         std::vector<double> &occupation);

  static void nonrelativistic_atomic_states(int core_charge,
                                            std::vector<int> &n,
                                            std::vector<int> &l,
                                            std::vector<double> &occupation);
};

}  // namespace lsms

#endif  // LSMS_STATES_HPP
