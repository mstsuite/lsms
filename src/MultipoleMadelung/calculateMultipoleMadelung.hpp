//
// Created by F.Moitzi on 19.01.2022.
//

#include "Main/SystemParameters.hpp"

#ifndef CALC_MULTIPOLE_MADELUNG_HPP
#define CALC_MULTIPOLE_MADELUNG_HPP

void calculateMultiMadelungMatrices(LSMSSystemParameters &lsms,
                                    CrystalParameters &crystal,
                                    LocalTypeInfo &local, int lmax = 0);

#endif