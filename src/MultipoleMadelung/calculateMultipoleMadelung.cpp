//
// Created by F.Moitzi on 19.01.2022.
//

#include "calculateMultipoleMadelung.hpp"

#include "Misc/Coeficients.hpp"
#include "MultipoleMadelung.hpp"

void calculateMultiMadelungMatrices(LSMSSystemParameters &lsms,
                                    CrystalParameters &crystal,
                                    LocalTypeInfo &local,
                                    int lmax) {
  if (lsms.global.iprint >= 0) {
    std::printf("Madelung calculations!\n");
  }

  lsms::MultipoleMadelung madelung(lsms, crystal, local, lmax);

}
