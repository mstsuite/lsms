//
// Created by F.Moitzi on 11.10.2022.
//

#include "setupCell.hpp"

#include "num_digits.hpp"

void lsms::setupCell(LSMSSystemParameters &lsms, CrystalParameters &crystal,
                     LocalTypeInfo &local) {
  crystal.omega = std::abs((crystal.bravais(1, 0) * crystal.bravais(2, 1) -
                            crystal.bravais(2, 0) * crystal.bravais(1, 1)) *
                               crystal.bravais(0, 2) +
                           (crystal.bravais(2, 0) * crystal.bravais(0, 1) -
                            crystal.bravais(0, 0) * crystal.bravais(2, 1)) *
                               crystal.bravais(1, 2) +
                           (crystal.bravais(0, 0) * crystal.bravais(1, 1) -
                            crystal.bravais(1, 0) * crystal.bravais(0, 1)) *
                               crystal.bravais(2, 2));

  double ws_radius =
      std::cbrt(crystal.omega * 3.0 / (4.0 * M_PI * crystal.num_atoms));
  double ws_volumes = crystal.omega / crystal.num_atoms;

  for (int i = 0; i < local.num_local; i++) {
    local.atom[i].rInscribed = ws_radius;
    local.atom[i].voronoi.rInscribedSphere = ws_radius;
    local.atom[i].voronoi.omegaInt = 0.0;

    if (lsms.fixRMT == 0) {
      local.atom[i].generateNewMesh = true;
    }

    local.atom[i].omegaMT = ws_volumes;
    local.atom[i].omegaWS = ws_volumes;
    local.atom[i].rws = ws_radius;
    local.atom[i].rmt = ws_radius;
  }
}

void lsms::printCell(LSMSSystemParameters &lsms, CrystalParameters &crystal,
                     LocalTypeInfo &local, LSMSCommunication &comm) {
  double total_volume = 0.0;
  for (int i = 0; i < local.num_local; i++) {
    total_volume += local.atom[i].omegaMT;
  }
  globalSum(comm, total_volume);

  int size = -1;
  size = std::max(size, num_digits(static_cast<int>(crystal.omega)));
  size = std::max(size, num_digits(static_cast<int>(total_volume)));

  if (lsms.global.iprint >= -1 && comm.rank == 0) {
    std::printf("%-12s = %*.4f [bohr^3]\n", "Cell volume", size, crystal.omega);
    std::printf("%-12s = %*.4f [bohr^3]\n", "ASA volume", size, total_volume);
  }
}
