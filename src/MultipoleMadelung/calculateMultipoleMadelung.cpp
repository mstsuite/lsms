//
// Created by F.Moitzi on 19.01.2022.
//

#include "calculateMultipoleMadelung.hpp"

#include "MultipoleMadelung.hpp"

void calculateMultiMadelungMatrices(LSMSSystemParameters &lsms,
                                    CrystalParameters &crystal,
                                    LocalTypeInfo &local) {

  std::vector<int> gindex(local.num_local);

  if (lsms.global.iprint >= 0) {
    std::printf("Multipole madelung calculations!\n");
  }

  lsms::MultipoleMadelung madelung(crystal.bravais, crystal.position,
                                   crystal.num_atoms, 0, local.global_id);

  for (auto j = 0; j < local.num_local; j++) {
    local.atom[j].madelungMatrix.resize(crystal.num_atoms);

    for (auto i = 0; i < crystal.num_atoms; i++) {
      local.atom[j].madelungMatrix[i] = madelung.getMadSum(i, j);
    }
  }
}
