//
// Created by F.Moitzi on 18.09.2021.
//

#include <gtest/gtest.h>

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <cmath>
#include <complex>
#include <vector>

#include "Complex.hpp"
#include "Real.hpp"
#include "associatedLegendreFunction.hpp"

namespace {

inline int plmIdx(int l, int m) { return l * (l + 1) / 2 + m; }

void spherical_harmonics(std::vector<Real> &vec, int lmax,
                         std::vector<Complex> &ylm, std::vector<Real> &plm) {
  using namespace std::complex_literals;

  auto q2 = vec[0] * vec[0] + vec[1] * vec[1];
  auto r = sqrt(q2 + vec[2] * vec[2]);
  auto q = sqrt(q2);

  Complex iphi = 1i * atan2(vec[1] / q, vec[0] / q);

  associatedLegendreFunctionNormalized(vec[2] / r, lmax, plm.data());

  int kl = 0;

  for (int l = 0; l <= lmax; l++) {
    ylm[kl] = plm[plmIdx(l, 0)];

    for (int m = 1; m <= l; m++) {
      ylm[kl + 1] = plm[plmIdx(l, m)] * std::exp(iphi * (Complex)m);
      ylm[kl - 1] = std::conj(ylm[kl + 1]) * std::pow(-1, m);
    }

    kl += (l + 1) * 2;
  }
}

TEST(SphericalHarmonics, SphHarmTest1) {
  int lmax = 8;

  std::vector<Real> plm((lmax + 1) * (lmax + 2) / 2, 0.0);
  std::vector<Complex> ylm((lmax + 1) * (lmax + 1), 0.0);

  std::vector<double> vec(3);

  vec[0] = -4.8830810748721376;
  vec[1] = -4.8830810748721376;
  vec[2] = -4.8830810748721376;

  spherical_harmonics(vec, lmax, ylm, plm);

  // (0.28209479177387831,0.0000000000000000)
  EXPECT_NEAR(0.28209479177387831, std::real(ylm[0]), 1.0e-10);
  EXPECT_NEAR(0.0, std::imag(ylm[0]), 1.0e-10);

  // (-0.19947114020071646,0.19947114020071652)
  EXPECT_NEAR(-0.19947114020071646, std::real(ylm[1]), 1.0e-10);
  EXPECT_NEAR(0.19947114020071652, std::imag(ylm[1]), 1.0e-10);

  // (-0.28209479177387831,-0.0000000000000000)
  EXPECT_NEAR(-0.28209479177387831, std::real(ylm[2]), 1.0e-10);
  EXPECT_NEAR(0.0, std::imag(ylm[2]), 1.0e-10);

  //
  //// (0.19947114020071646,0.19947114020071652)
  EXPECT_NEAR(0.19947114020071646, std::real(ylm[3]), 1.0e-10);
  EXPECT_NEAR(0.19947114020071652, std::imag(ylm[3]), 1.0e-10);
}
}  // namespace
