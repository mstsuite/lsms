//
// Created by F.Moitzi on 07.01.2023.
//

#ifndef LSMS_SRC_MAIN_UPDATECHARGEPOTENTIAL_HPP_
#define LSMS_SRC_MAIN_UPDATECHARGEPOTENTIAL_HPP_

#include "LSMSCommunication.hpp"
#include "Matrix.hpp"
#include "Real.hpp"
#include "SystemParameters.hpp"

namespace lsms {

void updateChargePotential(LSMSSystemParameters &lsms, LocalTypeInfo &local);

}  // namespace lsms

#endif  // LSMS_SRC_MAIN_UPDATECHARGEPOTENTIAL_HPP_
