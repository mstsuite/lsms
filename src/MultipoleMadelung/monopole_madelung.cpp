//
// Created by F.Moitzi on 06.01.2022.
//

#include "monopole_madelung.hpp"

#define _USE_MATH_DEFINES

#include <algorithm>
#include <cmath>

#include "integer_factors.hpp"
#include "lattice_utils.hpp"
#include "madelung_term.hpp"
#include "utils.hpp"

double lsms::scaling_factor(const lsms::matrix<double> &bravais, int lmax,
                            int max_iter, double fstep) {
  // Scaled bravais lattice
  lsms::matrix<double> r_brav(3, 3);
  // Scaled reciprocal lattice
  lsms::matrix<double> k_brav(3, 3);

  // Get the shortest axis
  auto a0 =
      std::sqrt(bravais(0, 0) * bravais(0, 0) + bravais(1, 0) * bravais(1, 0) +
                bravais(2, 0) * bravais(2, 0));

  auto a1 =
      std::sqrt(bravais(0, 1) * bravais(0, 1) + bravais(1, 1) * bravais(1, 1) +
                bravais(2, 1) * bravais(2, 1));

  auto a2 =
      std::sqrt(bravais(0, 2) * bravais(0, 2) + bravais(1, 2) * bravais(1, 2) +
                bravais(2, 2) * bravais(2, 2));

  double scaling_fac = std::min({a0, a1, a2});
  auto eta = 0.5 + 0.1 * std::max({a0, a1, a2}) / scaling_fac;
  scaling_fac /= 2.0 * M_PI;

  std::vector<int> nm(3);

  for (int i = 0; i <= max_iter; i++) {
    r_brav = bravais;
    r_brav.scale(1 / scaling_fac);
    k_brav = 0.0;

    reciprocal_lattice(r_brav, k_brav, scaling_fac);

    // Radius of real space truncation sphere
    nm = real_space_multiplication(r_brav, lmax, eta);
    auto rscut = rs_trunc_radius(r_brav, lmax, eta, nm);
#ifdef LSMS_DEBUG
    std::cout << nm[0] << " " << nm[1] << " " << nm[2] << std::endl;
#endif

    // Calculate number of lattice vectors
    auto nrslat = num_latt_vectors(r_brav, rscut, nm);

    // Radius of reciprocal space
    nm = reciprocal_space_multiplication(k_brav, lmax, eta);
    auto kncut = kn_trunc_radius(k_brav, lmax, eta, nm);
#ifdef LSMS_DEBUG
    std::cout << nm[0] << " " << nm[1] << " " << nm[2] << std::endl;
#endif

    // Calculate number of lattice vectors
    auto nknlat = num_latt_vectors(k_brav, kncut, nm);

#ifdef LSMS_DEBUG
    std::printf("ALAT: %f %f %f\n", r_brav(0, 0), r_brav(1, 1), r_brav(2, 2));
    std::printf("RS: %f KN: %f RSLAT: %d KNLAT: %d SC: %f\n", rscut, kncut,
                nrslat, nknlat, scaling_fac);
#endif

    if (nknlat < nrslat / 2) {
      scaling_fac = scaling_fac - fstep;
    } else if (nrslat < nknlat / 2) {
      scaling_fac = scaling_fac + fstep;
    } else {
      break;
    }
  }

  return scaling_fac;
}

int lsms::num_latt_vectors(const lsms::matrix<double> &brav, double cut,
                           const std::vector<int> &nm) {
  int number = 0;

  // To also include vectors on the boarder
  double vcut2 = cut * cut + 1e-6;

  std::vector<double> vn(3, 0.0);

  for (int x = -nm[0]; x <= nm[0]; x++) {
    for (int y = -nm[1]; y <= nm[1]; y++) {
      for (int z = -nm[2]; z <= nm[2]; z++) {
        vn[0] = x * brav(0, 0) + y * brav(0, 1) + z * brav(0, 2);
        vn[1] = x * brav(1, 0) + y * brav(1, 1) + z * brav(1, 2);
        vn[2] = x * brav(2, 0) + y * brav(2, 1) + z * brav(2, 2);

        auto norm = norm_sq(vn.begin(), vn.end());

        if (norm <= vcut2) {
          number++;
        }
      }
    }
  }

  return number;
}

std::vector<int> lsms::real_space_multiplication(
    const lsms::matrix<double> &brav, int lmax, double eta) {
  std::vector<int> nm(3);
  std::vector<double> r(3, 0.0);

  for (int i = 0; i < 3; i++) {
    r[i] = std::sqrt(brav(0, i) * brav(0, i) + brav(1, i) * brav(1, i) +
                     brav(2, i) * brav(2, i));
  }

  for (int i = 0; i < 3; i++) {
    nm[i] = 0;
    auto term = 1.0;
    while (term > 0.5 * lsms::EPSI) {
      nm[i]++;
      auto gamma = gamma_func(nm[i] * r[i] / eta, lmax);
      term = gamma[lmax] / std::pow((nm[i] * r[i] / 2.0), lmax + 1);
    }
  }

  return nm;
}

double lsms::rs_trunc_radius(const lsms::matrix<double> &brav, int lmax,
                             double eta, const std::vector<int> &nm) {
  auto cut = 0.0;

  std::vector<double> r(3, 0.0);

  for (int i = 0; i < 3; i++) {
    r[i] = sqrt(brav(0, i) * brav(0, i) + brav(1, i) * brav(1, i) +
                brav(2, i) * brav(2, i));
  }

  for (int i = -1; i <= 1; i++) {
    for (int idx = 0; idx < 3; idx++) {
      r[idx] = i * brav(idx, 0) * nm[0];
    }

    for (int j = -1; j <= 1; j++) {
      for (int idx = 0; idx < 3; idx++) {
        r[idx] = r[idx] + j * brav(idx, 1) * nm[1];
      }

      for (int k = -1; k <= 1; k++) {
        for (int idx = 0; idx < 3; idx++) {
          r[idx] = r[idx] + k * brav(idx, 2) * nm[2];
        }

        cut = std::max(cut, norm(r.begin(), r.end()));
      }
    }
  }

  return cut;
}

double lsms::kn_trunc_radius(const lsms::matrix<double> &brav, int lmax,
                             double eta, const std::vector<int> &nm) {
  auto cut = 0.0;

  std::vector<double> r(3, 0.0);

  for (int i = 0; i < 3; i++) {
    r[i] = sqrt(brav(0, i) * brav(0, i) + brav(1, i) * brav(1, i) +
                brav(2, i) * brav(2, i));
  }

  for (int i = -1; i <= 1; i++) {
    for (int idx = 0; idx < 3; idx++) {
      r[idx] = i * brav(idx, 0) * nm[0];
    }

    for (int j = -1; j <= 1; j++) {
      for (int idx = 0; idx < 3; idx++) {
        r[idx] = r[idx] + j * brav(idx, 1) * nm[1];
      }

      for (int k = -1; k <= 1; k++) {
        for (int idx = 0; idx < 3; idx++) {
          r[idx] = r[idx] + k * brav(idx, 2) * nm[2];
        }

        cut = std::max(cut, norm(r.begin(), r.end()));
      }
    }
  }

  return cut;
}

std::vector<int> lsms::reciprocal_space_multiplication(
    const lsms::matrix<double> &brav, int lmax, double eta) {
  std::vector<int> nm(3);
  std::vector<double> r(3, 0.0);

  for (int i = 0; i < 3; i++) {
    r[i] = brav(0, i) * brav(0, i) + brav(1, i) * brav(1, i) +
           brav(2, i) * brav(2, i);
  }

  auto fac = eta * eta / 4.0;

  for (int i = 0; i < 3; i++) {
    nm[i] = 0;
    auto term = 1.0;
    while (term > 0.5 * lsms::EPSI) {
      nm[i]++;
      auto rm = nm[i] * nm[i] * r[i];
      term = exp(-fac * rm) * std::pow(sqrt(rm), lmax - 2);
    }
  }

  return nm;
}

double lsms::calculate_eta(lsms::matrix<double> &bravais) {
  auto a0 = sqrt(bravais(0, 0) * bravais(0, 0) + bravais(1, 0) * bravais(1, 0) +
                 bravais(2, 0) * bravais(2, 0));

  auto a1 = sqrt(bravais(0, 1) * bravais(0, 1) + bravais(1, 1) * bravais(1, 1) +
                 bravais(2, 1) * bravais(2, 1));

  auto a2 = sqrt(bravais(0, 2) * bravais(0, 2) + bravais(1, 2) * bravais(1, 2) +
                 bravais(2, 2) * bravais(2, 2));

  auto scaling_fac = std::min({a0, a1, a2});

  return 0.5 + 0.1 * std::max({a0, a1, a2}) / scaling_fac;
}
