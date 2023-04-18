//
// Created by F.Moitzi on 19.01.2022.
//

#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"

#ifndef CALC_MULTIPOLE_MADELUNG_HPP
#define CALC_MULTIPOLE_MADELUNG_HPP

namespace lsms {

void calculateMultiMadelungMatrices(LSMSSystemParameters &lsms,
                                    CrystalParameters &crystal,
                                    LocalTypeInfo &local, int lmax = 0);


void printMultiMadelungMatrices(LSMSSystemParameters &lsms, LocalTypeInfo &local, LSMSCommunication &comm);


}

#endif