//
// Created by F.Moitzi on 03.01.2023.
//

#include "ChargeDensity.hpp"

#include "fmt/core.h"
#include "integrator.hpp"

namespace lsms {

void calculateRadialChargeDensity(LSMSSystemParameters &lsms,
                                  LocalTypeInfo &local) {
  /*
   * construct rho: rho(ir,0) = up density; rho(ir,1) = down density
   */

  for (int i = 0; i < local.num_local; i++) {
    auto &atom = local.atom[i];

    for (int is = 0; is < lsms.n_spin_pola; is++) {
      for (int ir = 0; ir < atom.jmt; ir++) {
        atom.rhoNew(ir, is) =
            -atom.r_mesh[ir] * atom.r_mesh[ir] * atom.greenint(ir, is) / M_PI;
      }

      int nMesh = atom.vr.l_dim();
      std::vector<Real> dens(nMesh);

      for (int ir = 0; ir < nMesh; ir++) {
        dens[ir] = atom.rhoNew(ir, is);
      }

      double norm = lsms::radialIntegral(dens, atom.r_mesh, atom.jmt);

      if (lsms.global.debug_radial_charge) {
        fmt::print("  int Qval [{}|{}]: {:40.24f}\n", i, is, norm);
      }
    }
  }

  /*
   * Add core and semi-core charge density
   */
  for (int i = 0; i < local.num_local; i++) {
    auto &atom = local.atom[i];

    for (int is = 0; is < lsms.n_spin_pola; is++) {
      for (int ir = 0; ir < atom.jmt; ir++) {
        atom.rhoNew(ir, is) += atom.corden(ir, is) + atom.semcor(ir, is);
      }

      if (lsms.global.debug_radial_charge) {
        int nMesh = atom.vr.l_dim();
        std::vector<Real> dens(nMesh);

        for (int ir = 0; ir < nMesh; ir++) {
          dens[ir] = atom.rhoNew(ir, is);
        }

        double norm = lsms::radialIntegral(dens, atom.r_mesh, atom.jmt);

        fmt::print("  int Qtot [{}|{}]: {:40.24f}\n", i, is, norm);
      }
    }
  }

}

void copyChargesAndPotential(LSMSSystemParameters &lsms,
                             LocalTypeInfo &local) {

  for (int i = 0; i < local.num_local; i++) {
    for (int is = 0; is < lsms.n_spin_pola; is++) {
      for (int ir = 0; ir < local.atom[i].r_mesh.size(); ir++) {
        local.atom[i].vr(ir, is) = local.atom[i].vrNew(ir, is);
        local.atom[i].rhotot(ir, is) = local.atom[i].rhoNew(ir, is);
      }
    }
  }

}

void checkRadialChargeDensity(LSMSSystemParameters &lsms,
                              LocalTypeInfo &local) {

  if (lsms.global.debug_radial_charge) {
    for (int i = 0; i < local.num_local; i++) {
      auto &atom = local.atom[i];

      for (int is = 0; is < lsms.n_spin_pola; is++) {

        int nMesh = atom.r_mesh.size();
        double norm = 0.0;
        std::vector<Real> dens(nMesh);

        for (int ir = 0; ir < nMesh; ir++) {
          dens[ir] = atom.rhoNew(ir, is);
        }

        norm = lsms::radialIntegral(dens, atom.r_mesh, atom.jmt);
        fmt::print("  int Qnew [{}|{}]: {:40.24f}\n", i, is, norm);

        for (int ir = 0; ir < nMesh; ir++) {
          dens[ir] = atom.rhotot(ir, is);
        }

        norm = lsms::radialIntegral(dens, atom.r_mesh, atom.jmt);
        fmt::print("  int Qold [{}|{}]: {:40.24f}\n", i, is, norm);

      }
    }

  }

}

void calculateCharge(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                     std::vector<Real> &qsub) {
  // Compute integrated densities of states and store in xval**
  // (from mufind_c.f)

  double total_charge = 0.0;

  for (int i = 0; i < local.num_local; i++) {
    if (lsms.n_spin_cant == 2)  // nspin == 3
    {
      local.atom[i].qvalmt = local.atom[i].dosckint[0];
      local.atom[i].qvalws = local.atom[i].dosint[0];
      local.atom[i].mvalws =
          local.atom[i].dosint[1] * local.atom[i].evecNew[0] +
              local.atom[i].dosint[2] * local.atom[i].evecNew[1] +
              local.atom[i].dosint[3] * local.atom[i].evecNew[2];
      local.atom[i].mvalmt =
          local.atom[i].dosckint[1] * local.atom[i].evecNew[0] +
              local.atom[i].dosckint[2] * local.atom[i].evecNew[1] +
              local.atom[i].dosckint[3] * local.atom[i].evecNew[2];
      local.atom[i].xvalmt[0] =
          0.5 * (local.atom[i].qvalmt + local.atom[i].mvalmt);
      local.atom[i].xvalwsNew[0] =
          0.5 * (local.atom[i].qvalws + local.atom[i].mvalws);
      local.atom[i].xvalmt[1] =
          0.5 * (local.atom[i].qvalmt - local.atom[i].mvalmt);
      local.atom[i].xvalwsNew[1] =
          0.5 * (local.atom[i].qvalws - local.atom[i].mvalws);
    } else if (lsms.n_spin_pola == 2)  // nspin = 2
    {
      local.atom[i].qvalmt =
          local.atom[i].dosckint[0] + local.atom[i].dosckint[1];
      local.atom[i].qvalws = local.atom[i].dosint[0] + local.atom[i].dosint[1];
      local.atom[i].mvalmt =
          local.atom[i].dosckint[0] - local.atom[i].dosckint[1];
      local.atom[i].mvalws = local.atom[i].dosint[0] - local.atom[i].dosint[1];

      local.atom[i].xvalmt[0] = local.atom[i].dosckint[0];
      local.atom[i].xvalwsNew[0] = local.atom[i].dosint[0];
      local.atom[i].xvalmt[1] = local.atom[i].dosckint[1];
      local.atom[i].xvalwsNew[1] = local.atom[i].dosint[1];
    } else  // nspin = 1
    {
      local.atom[i].qvalmt = local.atom[i].dosckint[0];
      local.atom[i].qvalws = local.atom[i].dosint[0];
      local.atom[i].mvalmt = 0.0;
      local.atom[i].mvalws = 0.0;

      local.atom[i].xvalmt[0] = local.atom[i].dosckint[0];
      local.atom[i].xvalwsNew[0] = local.atom[i].dosint[0];
    }

    local.atom[i].qtotws = 0.0;
    local.atom[i].qtotmt = 0.0;
    for (auto is = 0; is < lsms.n_spin_pola; is++) {
      local.atom[i].qtotws += local.atom[i].xvalwsNew[is] +
          local.atom[i].zsemss + local.atom[i].zcorss;
      local.atom[i].qtotmt += local.atom[i].xvalwsNew[is] +
          local.atom[i].zsemss + local.atom[i].zcorss;
    }

    local.atom[i].mtotws = local.atom[i].xvalwsNew[0] -
        local.atom[i].xvalwsNew[lsms.n_spin_pola - 1];

    if (lsms.global.debug_charge) {
      fmt::print("  exp Q [{}]:      {:40.24f}\n", i, local.atom[i].qtotws);
      total_charge += local.atom[i].qtotws;
    }

    qsub[local.global_id[i]] = local.atom[i].xvalwsNew[0] - local.atom[i].zvalss;
    if (lsms.n_spin_pola == 2) {
      qsub[local.global_id[i]] += local.atom[i].xvalwsNew[1];
    }

  }

  if (lsms.global.debug_charge) {
    fmt::print("  exp Qtotal:     {:40.24f}\n", total_charge);
  }
}

void calculateQRMS(LSMSSystemParameters &lsms, LocalTypeInfo &local) {
  for (int i = 0; i < local.num_local; i++) {
    std::vector<Real> rms(local.atom[i].jmt, 0.0);

    for (int is = 0; is < lsms.n_spin_pola; is++) {
      for (int ir = 0; ir < local.atom[i].jmt; ir++) {
        auto value =
            (local.atom[i].rhoNew(ir, is) - local.atom[i].rhotot(ir, is));
        rms[ir] = value * value;
      }

      double norm =
          lsms::radialIntegral(rms, local.atom[i].r_mesh, local.atom[i].jmt);

      local.atom[i].qrms[is] = std::sqrt(norm) / (local.atom[i].omegaMT);
    }
  }
}

}  // namespace lsms