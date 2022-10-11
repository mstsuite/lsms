//
// Created by F.Moitzi on 07.01.2023.
//

#include "updateChargePotential.hpp"

namespace lsms {

void updateChargePotential(LSMSSystemParameters &lsms, LocalTypeInfo &local) {
  for (int i = 0; i < local.num_local; i++) {
    local.atom[i].vdif = local.atom[i].vdifNew;

    for (int is = 0; is < lsms.n_spin_pola; is++) {
      for (int ir = 0; ir < local.atom[i].jmt; ir++) {
        local.atom[i].vr(ir, is) = local.atom[i].vrNew(ir, is);
        local.atom[i].rhotot(ir, is) = local.atom[i].rhoNew(ir, is);
      }

      local.atom[i].xvalwsNew[is] = local.atom[i].xvalws[is];
    }
  }
}

}  // namespace lsms