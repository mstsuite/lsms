//
// Created by F.Moitzi on 24.10.2022.
//

#include "newPotential.hpp"

#include <fmt/core.h>
#include <fmt/format.h>
#include <fmt/printf.h>

#include "Misc/integrator.hpp"
#include "Misc/poisson.hpp"

void lsms::calculateVRMS(LSMSSystemParameters &lsms, LocalTypeInfo &local) {
  for (int i = 0; i < local.num_local; i++) {
    std::vector<Real> rms(local.atom[i].jmt, 0.0);

    for (int is = 0; is < lsms.n_spin_pola; is++) {
      for (int ir = 0; ir < local.atom[i].jmt; ir++) {
        auto value = (local.atom[i].vrNew(ir, is) - local.atom[i].vr(ir, is));
        rms[ir] = value * value;
      }

      double norm =
          lsms::radialIntegral(rms, local.atom[i].r_mesh, local.atom[i].jmt);

      local.atom[i].vrms[is] = std::sqrt(norm) / (local.atom[i].omegaMT);
    }
  }
}

/**
 * @brief Calculate ASA potential
 * @param comm
 * @param lsms
 * @param local
 * @param crystal
 * @param qsub
 *
 * Code was cross check and potential generation should be fine
 *
 */
void lsms::ASAPotential::calculatePotential(LSMSCommunication &comm,
                                            LSMSSystemParameters &lsms,
                                            LocalTypeInfo &local,
                                            CrystalParameters &crystal,
                                            std::vector<Real> &qsub) {
  Real vpot_ave_hartree = 0.0;
  Real vpot_ave_core = 0.0;
  Real vpot_ave_xc = 0.0;
  Real vpot_ave_vmt1 = 0.0;

  Real vpot_ave = 0.0;
  Real u0Sum = 0.0;

  std::vector<Real> ro3(local.num_local, 0.0);
  std::vector<Real> dz(local.num_local, 0.0);

  for (int i = 0; i < local.num_local; i++) {
    /**
     * Site indices
     */

    int jmt;
    int jws;
    Real rSphere;

    jmt = local.atom[i].jmt;
    jws = local.atom[i].jws;
    rSphere = local.atom[i].rInscribed;

    /**
     * 1. Site-dependent and constant Madelung part of potential
     */

    Real vpot = 0.0;
    Real vmt1 = 0.0;
    Real u0 = 0.0;

    for (auto j = 0; j < crystal.num_atoms; j++) {

      int global_i = local.global_id[i];

      vmt1 -= local.atom[i].madelungMatrix[j] * qsub[j] * 2.0;

      u0 += local.atom[i].madelungMatrix[j] * qsub[j] * qsub[global_i];

    }

    local.atom[i].localMadelungEnergy = u0;
    local.atom[i].localMadelungPotential = vmt1;
    u0Sum += u0;

    /**
     * 2. Exchange-correlation potential
     */

    local.atom[i].vrNew = 0.0;
    local.atom[i].exchangeCorrelationV[0] = 0.0;  // vxcout
    local.atom[i].exchangeCorrelationV[1] = 0.0;
    local.atom[i].exchangeCorrelationE = 0.0;  // excout

    if (lsms.xcFunctional[0] == 0)  // built in functionals
    {
      int iexch = lsms.xcFunctional[1];  // built in functional

      for (auto is = 0; is < lsms.n_spin_pola; is++) {
        Real spin = 1.0 - Real(is) * 2.0;

        newexchg_(&lsms.n_spin_pola, &spin, &local.atom[i].rhoNew(0, 0),
                  &local.atom[i].rhoNew(0, lsms.n_spin_pola - 1),
                  &local.atom[i].exchangeCorrelationPotential(0, is),
                  &local.atom[i].exchangeCorrelationEnergy(0, is),
                  &local.atom[i].exchangeCorrelationV[is],
                  &local.atom[i].exchangeCorrelationE, &ro3[i], &dz[i],
                  local.atom[i].r_mesh.data(), &jmt, &iexch);
      }
    } else if (lsms.xcFunctional[0] == 1)  // libxc functionals
    {
#ifdef USE_LIBXC

      lsms.exch_corr->evaluate(local.atom[i].r_mesh, local.atom[i].h,
                               local.atom[i].rhoNew, jmt,
                               local.atom[i].exchangeCorrelationEnergy,
                               local.atom[i].exchangeCorrelationPotential);

#else
      fmt::printf("LSMS was not built with libxc support!!\n");
      MPI_Abort(comm.comm, 1);
#endif
    } else {
      fmt::printf("Unknown XC Functional class!\n");
      MPI_Abort(comm.comm, 1);
    }

    /**
     * 3. Hartree potential
     */

    std::vector<Real> vhartreederiv(local.atom[i].r_mesh.size(), 0.0);
    std::vector<Real> vhartree(local.atom[i].r_mesh.size(), 0.0);
    std::vector<Real> density(local.atom[i].r_mesh.size(), 0.0);

    if (lsms.n_spin_pola == 1) {
      for (auto ir = 0; ir < local.atom[i].r_mesh.size(); ir++) {
        density[ir] = local.atom[i].rhoNew(ir, 0);
      }
    } else {
      for (auto ir = 0; ir < local.atom[i].r_mesh.size(); ir++) {
        density[ir] =
            (local.atom[i].rhoNew(ir, 0) + local.atom[i].rhoNew(ir, 1));
      }
    }

    auto Vend = lsms::radialIntegral(density, local.atom[i].r_mesh, jmt) / rSphere;

    lsms::radial_poisson(vhartree, vhartreederiv, local.atom[i].r_mesh,
                         local.atom[i].h, density, jmt);

    /**
     * 4. Sum up all contributions
     */

    for (auto is = 0; is < lsms.n_spin_pola; is++) {
      for (auto ir = 0; ir < jmt; ir++) {
        vpot = 2.0 * vhartree[ir]
            - 2.0 * local.atom[i].ztotss / local.atom[i].r_mesh[ir]
            + local.atom[i].exchangeCorrelationPotential(ir, is)
            + vmt1;

        local.atom[i].vrNew(ir, is) = vpot * local.atom[i].r_mesh[ir];
      }

      vpot = 2.0 * vhartree[jmt - 1]
          - 2.0 * local.atom[i].ztotss / local.atom[i].r_mesh[jmt - 1]
          + local.atom[i].exchangeCorrelationPotential(jmt - 1, is)
          + vmt1;

      if (lsms.global.debug_potential) {

        std::stringstream ss;

        ss << fmt::sprintf(" ==== Potential debug (%d) (%d) ====\n", comm.rank, local.global_id[i]);
        ss << fmt::sprintf("  Vhartree: %30.24f  %30.24f\n", 2.0 * vhartree[jmt - 1], 2.0 * vhartree[0]);
        ss << fmt::sprintf("  Vcore:    %30.24f\n", -2.0 * local.atom[i].ztotss / local.atom[i].r_mesh[jmt - 1]);
        ss << fmt::sprintf("  VXC:      %30.24f  %30.24f\n", local.atom[i].exchangeCorrelationPotential(jmt - 1, is),
                           local.atom[i].exchangeCorrelationPotential(0, is));
        ss << fmt::sprintf("  VM:       %30.24f\n", vmt1);
        ss << fmt::sprintf("  Density   %30.20f  %30.24f\n",
                           local.atom[i].rhoNew(jmt - 1, 0),
                           local.atom[i].rhoNew(0, 0));
        ss << fmt::sprintf("  r         %30.20f  %30.24f\n", local.atom[i].r_mesh[jmt - 1], local.atom[i].r_mesh[0]);
        ss << fmt::sprintf("  VMT Vend: %30.24f  %30.24f\n", vmt1, 2.0 * Vend);

        std::cout << ss.str() << std::endl;
      }

      if (lsms.global.debug_potential) {

        vpot_ave_hartree += 2.0 * vhartree[jmt - 1];
        vpot_ave_core -= 2.0 * local.atom[i].ztotss / local.atom[i].r_mesh[jmt - 1];
        vpot_ave_xc += local.atom[i].exchangeCorrelationPotential(jmt - 1, is);
        vpot_ave_vmt1 += vmt1;

      }

      vpot_ave += vpot;
    }
  }

  /**
   * Calculate MTZ Zero
   */

  globalSum(comm, u0Sum);
  lsms.u0 = u0Sum;
  lsms.u0MT = u0Sum;

  globalSum(comm, vpot_ave);
  vpot_ave /= ((Real) lsms.n_spin_pola * lsms.num_atoms);

  if (lsms.global.debug_potential) {
    globalSum(comm, vpot_ave_hartree);
    globalSum(comm, vpot_ave_core);
    globalSum(comm, vpot_ave_xc);
    globalSum(comm, vpot_ave_vmt1);

    vpot_ave_hartree /= ((Real) lsms.n_spin_pola * lsms.num_atoms);
    vpot_ave_core /= ((Real) lsms.n_spin_pola * lsms.num_atoms);
    vpot_ave_xc /= ((Real) lsms.n_spin_pola * lsms.num_atoms);
    vpot_ave_vmt1 /= ((Real) lsms.n_spin_pola * lsms.num_atoms);
  }

  if (comm.rank == 0 && lsms.global.debug_potential) {

    std::stringstream ss;

    ss << fmt::sprintf(" ==== Potential debug (%d) ====\n", comm.rank);
    ss << fmt::sprintf("  Pot. ave.:         %40.24f\n", vpot_ave);
    ss << fmt::sprintf("  Pot. ave. Hartree: %40.24f\n", vpot_ave_hartree);
    ss << fmt::sprintf("  Pot. ave. Core:    %40.24f\n", vpot_ave_core);
    ss << fmt::sprintf("  Pot. ave. XC:      %40.24f\n", vpot_ave_xc);
    ss << fmt::sprintf("  Pot. ave. VMT1:    %40.24f\n", vpot_ave_vmt1);

    std::cout << ss.str() << std::endl;
  }

  for (auto i = 0; i < local.num_local; i++) {
    int jmt = local.atom[i].jmt;

    for (auto is = 0; is < lsms.n_spin_pola; is++) {
      for (auto ir = 0; ir < jmt; ir++) {
        local.atom[i].vrNew(ir, is) -= vpot_ave * local.atom[i].r_mesh[ir];
      }

      for (auto ir = jmt; ir < local.atom[i].r_mesh.size(); ir++) {
        local.atom[i].vrNew(ir, is) = 0.0;
      }
    }
  }

  lsms.vmt = vpot_ave;
}
