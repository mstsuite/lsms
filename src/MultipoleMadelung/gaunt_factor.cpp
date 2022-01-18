//
// Created by F.Moitzi on 04.01.2022.
//

#include "gaunt_factor.hpp"

double lsms::get_gaunt_factor(lsms::array3d<double> &cgnt,
                              lsms::matrix<int> &nj3, lsms::array3d<int> &kj3,
                              int k1, int k2, int k3) {
  double cg = 0.0;

  for (auto j3 = 0; j3 < nj3(k1, k2); j3++) {
    if (kj3(j3, k1, k2) == k3 + 1) {
      return cgnt(j3, k1, k2);
    }
  }

  return cg;
}
