//
// Created by F.Moitzi on 19.01.2022.
//

#include "calculateMultipoleMadelung.hpp"

#include "MultipoleMadelung.hpp"

#include "Misc/Coeficients.hpp"

void calculateMultiMadelungMatrices(LSMSSystemParameters &lsms,
                                    CrystalParameters &crystal,
                                    LocalTypeInfo &local) {
  if (lsms.global.iprint >= 0) {
    std::printf("Madelung calculations!\n");
  }

  lsms::MultipoleMadelung madelung(lsms, crystal, local);



}
