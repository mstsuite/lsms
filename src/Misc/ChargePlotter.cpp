//
// Created by F.Moitzi on 01.12.2022.
//

#include "ChargePlotter.hpp"

#include <numeric>

#include "integrator.hpp"
#include "mpi.h"

namespace lsms {

ChargePlotter::ChargePlotter(LSMSCommunication &comm,
                             LSMSSystemParameters &lsms,
                             std::string file_name) {
  if (comm.rank == 0) {
    file.open(file_name);
  }
};

void ChargePlotter::plotCharges(LSMSCommunication &comm,
                                LSMSSystemParameters &lsms,
                                LocalTypeInfo &local, int iter) {
  int local_num_atoms = local.num_local;
  int num_atoms = lsms.num_atoms;

  if (lsms.n_spin_pola == 2) {
    std::vector<double> qtot_up(local_num_atoms);
    std::vector<double> qtot_down(local_num_atoms);

    std::vector<double> global_qtot_up(num_atoms, 0.0);
    std::vector<double> global_qtot_down(num_atoms, 0.0);

    for (int i = 0; i < local_num_atoms; i++) {
      std::vector<double> density(local.atom[i].r_mesh.size());

      for (auto ir = 0; ir < local.atom[i].r_mesh.size(); ir++) {
        density[ir] = local.atom[i].rhoNew(ir, 0);
      }

      qtot_up[i] =
          radialIntegral(density, local.atom[i].r_mesh, local.atom[i].jmt);

      for (auto ir = 0; ir < local.atom[i].r_mesh.size(); ir++) {
        density[ir] = local.atom[i].rhoNew(ir, 1);
      }

      qtot_down[i] =
          radialIntegral(density, local.atom[i].r_mesh, local.atom[i].jmt);
    }

    // Global sum step

    for (int i = 0; i < local_num_atoms; i++) {
      int id = local.global_id[i];
      global_qtot_up[id] = qtot_up[i];
      global_qtot_down[id] = qtot_down[i];
    }

    if (comm.rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, global_qtot_down.data(), num_atoms, MPI_DOUBLE,
                 MPI_SUM, 0, comm.comm);
    } else {
      MPI_Reduce(global_qtot_down.data(), global_qtot_down.data(), num_atoms,
                 MPI_DOUBLE, MPI_SUM, 0, comm.comm);
    }

    if (comm.rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, global_qtot_up.data(), num_atoms, MPI_DOUBLE,
                 MPI_SUM, 0, comm.comm);
    } else {
      MPI_Reduce(global_qtot_up.data(), global_qtot_up.data(), num_atoms,
                 MPI_DOUBLE, MPI_SUM, 0, comm.comm);
    }

    auto result_up = std::reduce(global_qtot_up.begin(), global_qtot_up.end());

    // Print step
    file << std::setw(5) << iter << " ";
    for (int i = 0; i < num_atoms; i++) {
      file << std::fixed << std::showpoint << std::setprecision(12)
           << global_qtot_up[i] << " ";
    }
    file << std::fixed << std::showpoint << std::setprecision(12) << result_up
         << " ";

    auto result_down =
        std::reduce(global_qtot_down.begin(), global_qtot_down.end());

    for (int i = 0; i < num_atoms; i++) {
      file << std::fixed << std::showpoint << std::setprecision(12)
           << global_qtot_down[i] << " ";
    }
    file << std::fixed << std::showpoint << std::setprecision(12) << result_down
         << " ";
    file << std::endl;

  } else {
    std::vector<double> qtot_up(local_num_atoms);

    std::vector<double> global_qtot_up(num_atoms, 0.0);

    for (int i = 0; i < local_num_atoms; i++) {
      std::vector<double> density(local.atom[i].r_mesh.size());

      for (auto ir = 0; ir < local.atom[i].r_mesh.size(); ir++) {
        density[ir] = local.atom[i].rhoNew(ir, 0);
      }

      qtot_up[i] = local.atom[i].xvalwsNew[0];
      // radialIntegral(density, local.atom[i].r_mesh, local.atom[i].jmt);
    }

    // Global sum step

    for (int i = 0; i < local_num_atoms; i++) {
      int id = local.global_id[i];
      global_qtot_up[id] = qtot_up[i];
    }

    if (comm.rank == 0) {
      MPI_Reduce(MPI_IN_PLACE, global_qtot_up.data(), num_atoms, MPI_DOUBLE,
                 MPI_SUM, 0, comm.comm);
    } else {
      MPI_Reduce(global_qtot_up.data(), global_qtot_up.data(), num_atoms,
                 MPI_DOUBLE, MPI_SUM, 0, comm.comm);
    }

    auto result = std::reduce(global_qtot_up.begin(), global_qtot_up.end());

    // Print step
    file << std::setw(5) << iter << " ";
    for (int i = 0; i < num_atoms; i++) {
      file << std::fixed << std::showpoint << std::setprecision(12)
           << global_qtot_up[i] << " ";
    }
    file << result;
    file << std::endl;
  }
}

ChargePlotter::~ChargePlotter() {
  if (file.is_open()) {
    file.close();
  }
};

}  // namespace lsms