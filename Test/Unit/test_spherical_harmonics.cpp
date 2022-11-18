//
// Created by F.Moitzi on 18.09.2021.
//

#include <gtest/gtest.h>

#include <complex>
#include <vector>

#include "Complex.hpp"
#include "Real.hpp"
#include "SphericalHarmonics.hpp"
#include "accel_common.hpp"
#include "spherical_harmonics.hpp"

namespace spherical_harmonics_tests {

TEST(SphericalHarmonics, SphHarmTest2) {
  int lmax = 8;

  lsms::SphericalHarmonics sph(lmax);

  std::vector<Real> plm((lmax + 1) * (lmax + 2) / 2, 0.0);
  std::vector<Complex> ylm((lmax + 1) * (lmax + 1), 0.0);
  std::vector<Complex> ylm_compare((lmax + 1) * (lmax + 1), 0.0);

  std::vector<double> vec(3);

  vec[0] = 0.0;
  vec[1] = 0.0;
  vec[2] = 0.0;

  sph_harm_1_(vec.data(), &lmax, ylm.data());

  sph.computeYlm(lmax, vec, ylm_compare);

  EXPECT_NEAR(std::real(ylm[0]), std::real(ylm_compare[0]), 1.0e-12);
  EXPECT_NEAR(std::imag(ylm[0]), std::imag(ylm_compare[0]), 1.0e-12);

  EXPECT_NEAR(std::real(ylm[1]), std::real(ylm_compare[1]), 1.0e-12);
  EXPECT_NEAR(std::imag(ylm[1]), std::imag(ylm_compare[1]), 1.0e-12);
}

TEST(SphericalHarmonics, SphHarmTest3) {
  int lmax = 8;

  lsms::SphericalHarmonics sph(lmax);

  std::vector<Real> plm((lmax + 1) * (lmax + 2) / 2, 0.0);
  std::vector<Complex> ylm((lmax + 1) * (lmax + 1), 0.0);
  std::vector<Complex> ylm_compare((lmax + 1) * (lmax + 1), 0.0);

  std::vector<double> vec(3);

  vec[0] = -4.8830810748721376;
  vec[1] = -4.8830810748721376;
  vec[2] = -4.8830810748721376;

  sph_harm_1_(vec.data(), &lmax, ylm.data());

  sph.computeYlm(lmax, vec, ylm_compare);

  int k = 0;

  for (int l = 0; l <= lmax; l++) {
    for (int m = 1; m <= l; m++) {
      EXPECT_NEAR(std::real(ylm[k]), std::real(ylm_compare[k]), 1.0e-12);
      EXPECT_NEAR(std::imag(ylm[k]), std::imag(ylm_compare[k]), 1.0e-12);
      k++;
    }
  }
}

TEST(SphericalHarmonics, SphHarmTest4) {
  int lmax = 8;

  lsms::SphericalHarmonics sph(lmax);

  std::vector<Real> plm((lmax + 1) * (lmax + 2) / 2, 0.0);
  std::vector<Complex> ylm((lmax + 1) * (lmax + 1), 0.0);
  std::vector<Complex> ylm_compare((lmax + 1) * (lmax + 1), 0.0);

  std::vector<double> vec(3);

  vec[0] = 0.0;
  vec[1] = 0.0;
  vec[2] = -4.8830810748721376;

  sph_harm_1_(vec.data(), &lmax, ylm.data());

  sph.computeYlm(lmax, vec, ylm_compare);

  int k = 0;

  for (int l = 0; l <= lmax; l++) {
    for (int m = 1; m <= l; m++) {
      EXPECT_NEAR(std::real(ylm[k]), std::real(ylm_compare[k]), 1.0e-12);
      EXPECT_NEAR(std::imag(ylm[k]), std::imag(ylm_compare[k]), 1.0e-12);
      k++;
    }
  }
}

TEST(SphericalHarmonics, SphHarmTest5) {
  int lmax = 8;

  lsms::SphericalHarmonics sph(lmax);

  std::vector<Real> plm((lmax + 1) * (lmax + 2) / 2, 0.0);
  std::vector<Complex> ylm((lmax + 1) * (lmax + 1), 0.0);
  std::vector<Complex> ylm_compare((lmax + 1) * (lmax + 1), 0.0);

  std::vector<double> vec(3);

  vec[0] = 0.0;
  vec[1] = 0.0;
  vec[2] = 4.8830810748721376;

  sph_harm_1_(vec.data(), &lmax, ylm.data());

  sph.computeYlm(lmax, vec, ylm_compare);

  int k = 0;

  for (int l = 0; l <= lmax; l++) {
    for (int m = 1; m <= l; m++) {
      EXPECT_NEAR(std::real(ylm[k]), std::real(ylm_compare[k]), 1.0e-12);
      EXPECT_NEAR(std::imag(ylm[k]), std::imag(ylm_compare[k]), 1.0e-12);
      k++;
    }
  }
}

}  // namespace spherical_harmonics_tests
