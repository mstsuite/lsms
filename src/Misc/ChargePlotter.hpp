//
// Created by F.Moitzi on 01.12.2022.
//

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "LSMSCommunication.hpp"
#include "SystemParameters.hpp"

#ifndef LSMS_CHARGEPLOTTER_HPP
#define LSMS_CHARGEPLOTTER_HPP

namespace lsms {

class ChargePlotter {
  std::ofstream file;

 public:
  ChargePlotter(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                std::string file_name);

  ~ChargePlotter();

  void plotCharges(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                   LocalTypeInfo &local, int iter);
};

}  // namespace lsms

#endif  // LSMS_CHARGEPLOTTER_HPP
