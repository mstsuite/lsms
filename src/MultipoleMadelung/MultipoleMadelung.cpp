//
// Created by F.Moitzi on 16.12.2021.
//

#include "MultipoleMadelung.hpp"

#define _USE_MATH_DEFINES

#include <algorithm>
#include <complex>
#include <vector>

#include "integer_factors.hpp"
#include "lattice_utils.hpp"
#include "madelung.hpp"

lsms::MultipoleMadelung::MultipoleMadelung(
    const matrix<double> &lattice, const matrix<double> &atom_position,
    int num_atoms, int lmax, std::vector<int> global_position_index)
    : lmax{lmax}, num_atoms{num_atoms} {
  local_num_atoms =
      global_position_index
          .size();  // NOLINT(cppcoreguidelines-narrowing-conversions)

  int kmax = get_kmax(lmax);
  int jmax = get_jmax(lmax);

  auto r_brav = lattice;
  auto k_brav = lattice;

  // 1. Scaling factors and rescale lattices and atomic positions
  scaling_factor = lsms::scaling_factor(lattice, lmax);
  r_brav.scale(1.0 / scaling_factor);

  // atom_position.scale(1.0 / scaling_factor);

  reciprocal_lattice(r_brav, k_brav, scaling_factor);

  // auto omegbra = lsms::omega(r_brav);

  // 2. Calculate truncation spheres
  auto eta = lsms::calculate_eta(r_brav);

  // Real-space
  r_nm = lsms::real_space_multiplication(r_brav, lmax, eta);
  rscut = lsms::rs_trunc_radius(r_brav, lmax, eta, r_nm);
  nrslat = num_latt_vectors(r_brav, rscut, r_nm);

  // Reciprocal-space
  k_nm = lsms::reciprocal_space_multiplication(k_brav, lmax, eta);
  kncut = lsms::kn_trunc_radius(k_brav, lmax, eta, k_nm);
  nknlat = num_latt_vectors(k_brav, kncut, k_nm);

#ifdef LSMS_DEBUG
  std::printf("%d %d %d: %lf %d\n", r_nm[0], r_nm[1], r_nm[2], rscut, nrslat);
  std::printf("%d %d %d: %lf %d\n", k_nm[0], k_nm[1], k_nm[2], kncut, nknlat);
#endif

  // 3. Create the lattices
  matrix<double> rslat;
  std::vector<double> rslatsq;

  matrix<double> knlat;
  std::vector<double> knlatsq;

  std::tie(rslat, rslatsq) =
      lsms::create_lattice_and_sq(r_brav, rscut, r_nm, nrslat);

  std::tie(knlat, knlatsq) =
      lsms::create_lattice_and_sq(k_brav, kncut, k_nm, nknlat);

  // 4. Calculate the Madelung matrix and the prefactor matrix
  madsum = matrix<double>(num_atoms, local_num_atoms);

  if (jmax > 1) {
    dl_matrix = array3d<std::complex<double>>(num_atoms, kmax, local_num_atoms);
  }

  for (auto i{0}; i < global_position_index.size(); i++) {
    lsms::calculate_madelung_matrix(num_atoms, global_position_index[i], i,
                                    lmax, eta, scaling_factor, r_brav,
                                    atom_position, rslat, rslatsq, knlat,
                                    knlatsq, madsum, dl_matrix);
  }

  // 5. Dl factors
  if (jmax > 1) {
    dl_factor = lsms::calculate_dl_factor(lmax);
  }
}

double lsms::MultipoleMadelung::getMadSum(int i, int j) const {
  return madsum(i, j);
}

std::complex<double> lsms::MultipoleMadelung::getDlMatrix(int i, int k, int j) {
  return dl_matrix(i, k, j);
}

double lsms::MultipoleMadelung::getDlFactor(int i, int j) const {
  return dl_factor(i, j);
}

double lsms::MultipoleMadelung::getScalingFactor() const {
  return scaling_factor;
}

double lsms::MultipoleMadelung::getRsCut() const { return rscut; }

double lsms::MultipoleMadelung::getKnCut() const { return kncut; }

std::vector<int> lsms::MultipoleMadelung::getKnSize() const { return k_nm; }

std::vector<int> lsms::MultipoleMadelung::getRsSize() const { return r_nm; }
