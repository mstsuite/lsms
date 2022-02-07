//
// Created by F.Moitzi on 25.01.2022.
//

#include "GauntFactor.hpp"

#include "integer_factors.hpp"

#include "gaunt_factor.hpp"

lsms::math::GauntFactor::GauntFactor(int lmax) {

  std::size_t kmax_mad = lsms::get_kmax(lmax);

  cgnt = lsms::NDArray<double, 3>(lmax + 1, (lmax + 1) * (lmax + 1),
                                  (lmax + 1) * (lmax + 1));

  nj3 = lsms::NDArray<int, 2>(kmax_mad, kmax_mad);

  kj3 = lsms::NDArray<int, 3>(lmax + 1, kmax_mad, kmax_mad);




  gaunt_factor(&lmax, cgnt.data(), kj3.data(), nj3.data());


}

double lsms::math::GauntFactor::table(int k1, int k2, int k3) const {

  double cg = 0.0;

  for (auto j3 = 0; j3 < nj3(k1, k2); j3++) {
    if (kj3(j3, k1, k2) == k3 + 1) {
      return cgnt(j3, k1, k2);
    }
  }

  return cg;
}
