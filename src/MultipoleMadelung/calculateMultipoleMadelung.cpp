//
// Created by F.Moitzi on 19.01.2022.
//

#include "calculateMultipoleMadelung.hpp"

#include "fmt/core.h"
#include "fmt/ranges.h"

#include "Misc/Coeficients.hpp"
#include "MultipoleMadelung.hpp"

void lsms::calculateMultiMadelungMatrices(LSMSSystemParameters &lsms,
                                          CrystalParameters &crystal,
                                          LocalTypeInfo &local, int lmax) {
  if (lsms.global.iprint >= 0) {
    std::printf("Madelung calculations!\n");
  }

  lsms::MultipoleMadelung madelung(lsms, crystal, local, lmax);
}

void lsms::printMultiMadelungMatrices(LSMSSystemParameters &lsms, LocalTypeInfo &local, LSMSCommunication &comm) {

  std::stringstream ss;

  ss << " Rank " << comm.rank << ": " << std::endl;

  for (auto i = 0; i < local.num_local; i++) {

    ss << std::endl;

    for (auto j = 0; j < lsms.num_atoms; j++) {
      ss << fmt::format("    {:3d} - {:3d}: {:24.16f}\n", local.global_id[i], j, local.atom[i].madelungMatrix[j] * 2.0);
    }

    ss << std::endl;

  }

  std::cout << ss.str() << std::endl;

}
