//
// Created by F.Moitzi on 21.12.2021.
//

#include "lattice_utils.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "utils.hpp"

double lsms::omega(lsms::matrix<double> &bravais) {
  return std::abs(
      (bravais(1, 0) * bravais(2, 1) - bravais(2, 0) * bravais(1, 1)) *
          bravais(0, 2) +
      (bravais(2, 0) * bravais(0, 1) - bravais(0, 0) * bravais(2, 1)) *
          bravais(1, 2) +
      (bravais(0, 0) * bravais(1, 1) - bravais(1, 0) * bravais(0, 1)) *
          bravais(2, 2));
}

void lsms::insert_ordered(lsms::matrix<double> &latt_vec,
                          std::vector<double> &latt_vec_sq, int len,
                          std::vector<double> &vec, double &v_sq) {
  for (int i = 0; i < len; i++) {
    if (latt_vec_sq[i] >= v_sq) {
      for (int j = len - 1; j >= i; j--) {
        latt_vec(0, j + 1) = latt_vec(0, j);
        latt_vec(1, j + 1) = latt_vec(1, j);
        latt_vec(2, j + 1) = latt_vec(2, j);
        latt_vec_sq[j + 1] = latt_vec_sq[j];
      }

      latt_vec_sq[i] = v_sq;

      latt_vec(0, i) = vec[0];
      latt_vec(1, i) = vec[1];
      latt_vec(2, i) = vec[2];

      return;
    }
  }

  latt_vec_sq[len] = v_sq;
  latt_vec(0, len) = vec[0];
  latt_vec(1, len) = vec[1];
  latt_vec(2, len) = vec[2];
}

__attribute__((unused)) lsms::matrix<double> lsms::create_lattice(
    const lsms::matrix<double> &brav, double cutoff, const std::vector<int> &nm,
    int size) {
  lsms::matrix<double> latt_vec(3, size);
  std::vector<double> latt_vec_sq(size);

  int number = 0;

  // To also include vectors on the boarder
  double vcut2 = cutoff * cutoff + 1e-6;

  std::vector<double> vn(3, 0.0);

  for (int x = -nm[0]; x <= nm[0]; x++) {
    for (int y = -nm[1]; y <= nm[1]; y++) {
      for (int z = -nm[2]; z <= nm[2]; z++) {
        vn[0] = x * brav(0, 0) + y * brav(0, 1) + z * brav(0, 2);
        vn[1] = x * brav(1, 0) + y * brav(1, 1) + z * brav(1, 2);
        vn[2] = x * brav(2, 0) + y * brav(2, 1) + z * brav(2, 2);

        auto norm = norm_sq(vn.begin(), vn.end());

        if (norm <= vcut2) {
          insert_ordered(latt_vec, latt_vec_sq, number, vn, norm);
          number++;
        }
      }
    }
  }

  return latt_vec;
}

std::tuple<lsms::matrix<double>, std::vector<double>>
lsms::create_lattice_and_sq(lsms::matrix<double> &brav, double cutoff,
                            const std::vector<int> &nm, int size) {
  lsms::matrix<double> latt_vec(3, size);
  std::vector<double> latt_vec_sq(size);

  int number = 0;

  // To also include vectors on the boarder
  double vcut2 = cutoff * cutoff + 1e-6;

  std::vector<double> vn(3, 0.0);

  for (int x = -nm[0]; x <= nm[0]; x++) {
    for (int y = -nm[1]; y <= nm[1]; y++) {
      for (int z = -nm[2]; z <= nm[2]; z++) {
        vn[0] = x * brav(0, 0) + y * brav(0, 1) + z * brav(0, 2);
        vn[1] = x * brav(1, 0) + y * brav(1, 1) + z * brav(1, 2);
        vn[2] = x * brav(2, 0) + y * brav(2, 1) + z * brav(2, 2);

        auto norm = norm_sq(vn.begin(), vn.end());

        if (norm <= vcut2) {
          insert_ordered(latt_vec, latt_vec_sq, number, vn, norm);
          number++;
        }
      }
    }
  }

  return std::make_tuple(latt_vec, latt_vec_sq);
}

void lsms::reciprocal_lattice(lsms::matrix<double> &bravais,
                              lsms::matrix<double> &reciprocal_bravais,
                              double &scale) {
  auto vol = (bravais(1, 0) * bravais(2, 1) - bravais(2, 0) * bravais(1, 1)) *
                 bravais(0, 2) +
             (bravais(2, 0) * bravais(0, 1) - bravais(0, 0) * bravais(2, 1)) *
                 bravais(1, 2) +
             (bravais(0, 0) * bravais(1, 1) - bravais(1, 0) * bravais(0, 1)) *
                 bravais(2, 2);

  auto factor = 2.0 * M_PI / fabs(vol);

  reciprocal_bravais(0, 0) =
      factor * (bravais(1, 1) * bravais(2, 2) - bravais(2, 1) * bravais(1, 2));
  reciprocal_bravais(1, 0) =
      factor * (bravais(2, 1) * bravais(0, 2) - bravais(0, 1) * bravais(2, 2));
  reciprocal_bravais(2, 0) =
      factor * (bravais(0, 1) * bravais(1, 2) - bravais(1, 1) * bravais(0, 2));

  reciprocal_bravais(0, 1) =
      factor * (bravais(1, 2) * bravais(2, 0) - bravais(2, 2) * bravais(1, 0));
  reciprocal_bravais(1, 1) =
      factor * (bravais(2, 2) * bravais(0, 0) - bravais(0, 2) * bravais(2, 0));
  reciprocal_bravais(2, 1) =
      factor * (bravais(0, 2) * bravais(1, 0) - bravais(1, 2) * bravais(0, 0));

  reciprocal_bravais(0, 2) =
      factor * (bravais(1, 0) * bravais(2, 1) - bravais(2, 0) * bravais(1, 1));
  reciprocal_bravais(1, 2) =
      factor * (bravais(2, 0) * bravais(0, 1) - bravais(0, 0) * bravais(2, 1));
  reciprocal_bravais(2, 2) =
      factor * (bravais(0, 0) * bravais(1, 1) - bravais(1, 0) * bravais(0, 1));
}