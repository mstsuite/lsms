//
// Created by F.Moitzi on 25.01.2022.
//

#include "GauntFactor.hpp"

#include "gaunt_factor.hpp"
#include "integer_factors.hpp"

lsms::math::GauntFactor::GauntFactor(int lmax) {
  std::size_t kmax_mad = lsms::get_kmax(lmax);

  cgnt = lsms::NDArray<double, 3>(lmax + 1, (lmax + 1) * (lmax + 1),
                                  (lmax + 1) * (lmax + 1));

  nj3 = lsms::NDArray<int, 2>(kmax_mad, kmax_mad);

  kj3 = lsms::NDArray<int, 3>(lmax + 1, kmax_mad, kmax_mad);

  gaunt_factor(&lmax, cgnt.data(), kj3.data(), nj3.data());
}

double lsms::math::GauntFactor::table(int l1, int m1, int l2, int m2, int l3,
                                      int m3) const {
  double cg = 0.0;

  auto k1 = l1 * l1 + l1 - m1;

  auto k3 = l2 * l2 + l2 - m2;

  auto k2 = l3 * l3 + l3 - m3;

  for (auto j3 = 0; j3 < nj3(k1, k2); j3++) {
    if (kj3(j3, k1, k2) == k3 + 1) {
      return cgnt(j3, k1, k2);
    }
  }


  return cg;
}

double lsms::math::GauntFactorWrapper::table(int l1, int m1,
                                             int l2, int m2,
                                             int l3, int m3) const {
  double cg = 0.0;

  auto k2 = l2 * l2 + l2 - m2;
  auto k3 = l3 * l3 + l3 - m3;

  auto w1 = l1 / 2;

  if (m1 + m2 == m3) {
    cg = cgnt(w1, k2, k3);
  }


  return cg;
}