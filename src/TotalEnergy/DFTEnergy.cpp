//
// Created by F.Moitzi on 05.07.2022.
//

#include "DFTEnergy.hpp"

#include <mpi.h>

template<class T>
static int num_digits(T number) {
  int digits = 0;

  if (number < 0) digits = 1; // remove this line if '-' counts as a digit

  while (number) {
    number /= 10;
    digits++;
  }
  return digits;
}

void lsms::print_dft_energy(const DFTEnergy &energy) {

  int size = -1;

  size = std::max(size, num_digits(static_cast<int> (energy.core_eigen)));
  size = std::max(size, num_digits(static_cast<int> (energy.semicore_eigen)));
  size = std::max(size, num_digits(static_cast<int> (energy.one_ele)));
  size = std::max(size, num_digits(static_cast<int> (energy.ks)));
  size = std::max(size, num_digits(static_cast<int> (energy.kinetic)));
  size = std::max(size, num_digits(static_cast<int> (energy.hartree)));
  size = std::max(size, num_digits(static_cast<int> (energy.core_hartree)));
  size = std::max(size, num_digits(static_cast<int> (energy.coloumb)));
  size = std::max(size, num_digits(static_cast<int> (energy.kinetic)));
  size = std::max(size, num_digits(static_cast<int> (energy.hartree)));
  size = std::max(size, num_digits(static_cast<int> (energy.core_hartree)));
  size = std::max(size, num_digits(static_cast<int> (energy.coloumb)));
  size = std::max(size, num_digits(static_cast<int> (energy.xc)));
  size = std::max(size, num_digits(static_cast<int> (energy.zero_point)));
  size = std::max(size, num_digits(static_cast<int> (energy.lsf)));
  size = std::max(size, num_digits(static_cast<int> (energy.madelung)));
  size = std::max(size, num_digits(static_cast<int> (energy.it_madelung)));
  size = std::max(size, num_digits(static_cast<int> (energy.it_xc)));
  size = std::max(size, num_digits(static_cast<int> (energy.total)));

  size += 12;

  std::printf("===================\n");

  std::printf("%-12s = %*.10f Ry\n", "Core", size, energy.core_eigen);
  std::printf("%-12s = %*.10f Ry\n", "Semicore", size, energy.semicore_eigen);
  std::printf("%-12s = %*.10f Ry\n", "One electron", size, energy.one_ele);
  std::printf("%-12s = %*.10f Ry\n", "Kohn-Sham", size, energy.ks);
  std::printf("%-12s = %*.10f Ry\n", "Kinetic", size, energy.kinetic);

  std::printf("%-12s = %*.10f Ry\n", "Hartree", size, energy.hartree);
  std::printf("%-12s = %*.10f Ry\n", "Core Hartree", size, energy.core_hartree);
  std::printf("%-12s = %*.10f Ry\n", "Coloumb", size, energy.coloumb);

  std::printf("%-12s = %*.10f Ry\n", "XC", size, energy.xc);
  std::printf("%-12s = %*.10f Ry\n", "ZPE", size, energy.zero_point);
  std::printf("%-12s = %*.10f Ry\n", "LSF", size, energy.lsf);

  std::printf("%-12s = %*.10f Ry\n", "MT Madelung", size, energy.madelung);

  std::printf("%-12s = %*.10f Ry\n", "IT Madelung", size, energy.it_madelung);
  std::printf("%-12s = %*.10f Ry\n", "IT XC", size, energy.it_xc);

  std::printf("%-12s = %*.10f Ry\n", "Total energy", size, energy.total);

  std::printf("===================\n");

}

static void sum_dft_energies(lsms::DFTEnergy *in, lsms::DFTEnergy *inout,
                             int *len, MPI_Datatype *type) {
  for (int i = 0; i < *len; i++) {
    inout[i] += in[i];
  }
}

void lsms::globalSum(LSMSCommunication &comm, DFTEnergy &dft_energy) {
  DFTEnergy global_dft_energy;

  MPI_Datatype dft_type;

  constexpr int size = 15;

  std::vector<MPI_Datatype> types(size, MPI_DOUBLE);
  std::vector<int> blocks(size, 1);
  std::vector<MPI_Aint> displacements(size);
  MPI_Aint base;

  MPI_Get_address(&dft_energy, &base);
  MPI_Get_address(&dft_energy.zero_point, &displacements[0]);
  MPI_Get_address(&dft_energy.core_eigen, &displacements[1]);
  MPI_Get_address(&dft_energy.semicore_eigen, &displacements[2]);
  MPI_Get_address(&dft_energy.one_ele, &displacements[3]);
  MPI_Get_address(&dft_energy.ks, &displacements[4]);
  MPI_Get_address(&dft_energy.kinetic, &displacements[5]);
  MPI_Get_address(&dft_energy.hartree, &displacements[6]);
  MPI_Get_address(&dft_energy.core_hartree, &displacements[7]);
  MPI_Get_address(&dft_energy.coloumb, &displacements[8]);
  MPI_Get_address(&dft_energy.xc, &displacements[9]);
  MPI_Get_address(&dft_energy.lsf, &displacements[10]);
  MPI_Get_address(&dft_energy.total, &displacements[11]);
  MPI_Get_address(&dft_energy.madelung, &displacements[12]);
  MPI_Get_address(&dft_energy.it_madelung, &displacements[13]);
  MPI_Get_address(&dft_energy.it_xc, &displacements[14]);

  for (auto &value: displacements) {
    value -= base;
  }

  MPI_Type_create_struct(size, blocks.data(), displacements.data(),
                         types.data(), &dft_type);

  MPI_Type_commit(&dft_type);

  MPI_Op dft_energy_sum;
  MPI_Op_create(
      reinterpret_cast<void (*)(void *, void *, int *, MPI_Datatype *)>(
          sum_dft_energies),
      1, &dft_energy_sum);

  MPI_Allreduce(&dft_energy, &global_dft_energy, 1, dft_type, dft_energy_sum,
                comm.comm);

  MPI_Op_free(&dft_energy_sum);

  dft_energy = global_dft_energy;

  MPI_Type_free(&dft_type);
}