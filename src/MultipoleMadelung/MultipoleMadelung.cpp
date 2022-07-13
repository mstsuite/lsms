//
// Created by F.Moitzi on 16.12.2021.
//

#include "MultipoleMadelung.hpp"

#include <complex>
#include <vector>

#include "lattice_utils.hpp"
#include "madelung_term.hpp"
#include "monopole_madelung.hpp"

lsms::MultipoleMadelung::MultipoleMadelung(LSMSSystemParameters &lsms,
                                           CrystalParameters &crystal,
                                           LocalTypeInfo &local)
    : num_atoms{crystal.num_atoms}, lmax{0}, local_num_atoms{local.num_local} {
  auto r_brav = crystal.bravais;
  auto k_brav = crystal.bravais;

  // Scaling factors for the optimal truncation sphere
  scaling_factor = lsms::scaling_factor(crystal.bravais, lmax);

  r_brav.scale(1.0 / scaling_factor);
  reciprocal_lattice(r_brav, k_brav, scaling_factor);

  // Calculate truncation spheres
  auto eta = lsms::calculate_eta(r_brav);

  // Real-space cutoffs
  double timeRealSpace = MPI_Wtime();
  r_nm = lsms::real_space_multiplication(r_brav, lmax, eta);
  rscut = lsms::rs_trunc_radius(r_brav, lmax, eta, r_nm);
  nrslat = num_latt_vectors(r_brav, rscut, r_nm);
  timeRealSpace = MPI_Wtime() - timeRealSpace;

  // Reciprocal-space cutoffs
  double timeReciprocalSpace = MPI_Wtime();
  k_nm = lsms::reciprocal_space_multiplication(k_brav, lmax, eta);
  kncut = lsms::kn_trunc_radius(k_brav, lmax, eta, k_nm);
  nknlat = num_latt_vectors(k_brav, kncut, k_nm);
  timeReciprocalSpace = MPI_Wtime() - timeReciprocalSpace;

  if (lsms.global.iprint >= 0) {
    std::printf("Eta: %lf Scaling: %lf\n", eta, scaling_factor);
    std::printf("Real space: %3d %3d %3d: %lf %8d\n", r_nm[0], r_nm[1], r_nm[2],
                rscut, nrslat);
    std::printf("Reciprocal space: %3d %3d %3d: %lf %8d\n", k_nm[0], k_nm[1],
                k_nm[2], kncut, nknlat);
  }

  // Create the lattice vectors
  matrix<double> rslat;
  std::vector<double> rslatsq;

  matrix<double> knlat;
  std::vector<double> knlatsq;

  std::tie(rslat, rslatsq) =
      lsms::create_lattice_and_sq(r_brav, rscut, r_nm, nrslat);

  std::tie(knlat, knlatsq) =
      lsms::create_lattice_and_sq(k_brav, kncut, k_nm, nknlat);

  /*
   * Calculate Madelung matrix for monopoles (lmax = 0)
   */
  auto omega = lsms::omega(r_brav);


  // Zero terms
  auto term0 = -M_PI * eta * eta / omega;

  for (auto j = 0; j < local.num_local; j++) {
    local.atom[j].madelungMatrix.resize(crystal.num_atoms);
  }

  double timeLoopSpace = MPI_Wtime();

  // Introduced for smaller objects
  auto position = crystal.position;

  #pragma omp parallel for collapse(2) firstprivate(nknlat, nrslat, scaling_factor, position, knlatsq, knlat, rslat, rslatsq) default(shared)
  for (int atom_i = 0; atom_i < num_atoms; atom_i++) {
    for (int local_i = 0; local_i < local_num_atoms; local_i++) {

      std::vector<double> aij(3);
      double r0tm;
      int ibegin;

      // Global index
      int global_i = local.global_id[local_i];

      // a_ij in unit of a0
      for (int idx = 0; idx < 3; idx++) {
        aij[idx] = position(idx, atom_i) / scaling_factor -
            position(idx, global_i) / scaling_factor;
      }

      // Real space terms: first terms
      if (global_i == atom_i) {
        ibegin = 1;
        r0tm = -2.0 / std::sqrt(M_PI) / eta;
      } else {
        ibegin = 0;
        r0tm = 0.0;
      }

      // Reciprocal space term
      auto term1 =
          reciprocal_space_term(knlat, knlatsq, aij, nknlat, eta, omega);

      // Real space term
      auto term2 = real_space_term(rslat, aij, nrslat, ibegin, eta);

      // Madelung matrix contribution
      local.atom[local_i].madelungMatrix[atom_i] =
          (term1 + term2 + r0tm + term0) / scaling_factor;
    }
  }

  timeLoopSpace = MPI_Wtime() - timeLoopSpace;

  if (lsms.global.iprint >= 0) {
    std::printf("Time: %16s %lf\n", "Real:", timeRealSpace);
    std::printf("Time: %16s %lf\n", "Reciprocal:", timeReciprocalSpace);
    std::printf("Time: %16s %lf\n", "Loop:", timeLoopSpace);
  }

}

double lsms::MultipoleMadelung::getScalingFactor() const {
  return scaling_factor;
}

double lsms::MultipoleMadelung::getRsCut() const { return rscut; }

double lsms::MultipoleMadelung::getKnCut() const { return kncut; }

std::vector<int> lsms::MultipoleMadelung::getKnSize() const { return k_nm; }

std::vector<int> lsms::MultipoleMadelung::getRsSize() const { return r_nm; }
