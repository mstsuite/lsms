//
// Created by F.Moitzi on 18.09.2021.
//

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "Coeficients.hpp"
#include "Misc/integrateOneDim.hpp"
#include "common.hpp"
#include "integrator.hpp"

using namespace lsms;

template <typename T>
bool approx_equal(T x, T y, T epsilon) {
  return fabs(x - y) / max(fabs(x), fabs(y)) <= epsilon;
}

template <typename T>
T relative_diff(T ref, T val) {
  return std::fabs(ref - val) / std::fabs(ref);
}

extern "C" {
void cgaunt_(int *lmax, Real *clm, Real *plmg, Real *tg, Real *wg, Real *cgnt,
             int *lofk, int *mofk, int *iprint, char *istop);
}

TEST(GauntFactorTest, Lmax3) { /*
                                * Reference results:
                                *
                                * MaxJ3 = lmax + 1
                                *
                                * cgnt(1:MaxJ3, 1:kmax, 1:kmax)
                                *
                                * cgnt(1, 1, 1)
                                * cgnt(1, 2, 1)
                                * cgnt(1, 1, 2)
                                * cgnt(1, 2, 2)
                                * cgnt(2, 2, 2)
                                * cgnt(3, 3, 3)
                                *
                                *   0.28209479177387842
                                *  -0.28209479177387836
                                *   0.28209479177387836
                                *  -0.12615662610100786
                                *   0.28209479177387836
                                *   0.0000000000000000
                                *
                                */

  int lmax = 3;

  int iprint = 2;
  char istop[32];

  AngularMomentumIndices a;
  a.init(lmax * 2);

  SphericalHarmonicsCoeficients s;
  s.init(lmax * 2);

  std::vector<Real> tg, wg;
  Matrix<Real> plmg;
  Array3d<double> cgnt;

  tg.resize(2 * (2 * lmax + 1));
  wg.resize(2 * (2 * lmax + 1));
  cgnt.resize(lmax + 1, (lmax + 1) * (lmax + 1), (lmax + 1) * (lmax + 1));
  plmg.resize(((2 * lmax + 1) * (2 * lmax + 2)) / 2, 2 * lmax + 1);

  cgaunt_(&lmax, &s.clm[0], &plmg(0, 0), &tg[0], &wg[0], &cgnt(0, 0, 0),
          &a.lofk[0], &a.mofk[0], &iprint, &istop[0]);

  EXPECT_NEAR(0.28209479177387842, cgnt(0, 0, 0), 1e-12);
  EXPECT_NEAR(-0.28209479177387842, cgnt(0, 1, 0), 1e-12);
  EXPECT_NEAR(0.28209479177387842, cgnt(0, 0, 1), 1e-12);
  EXPECT_NEAR(-0.12615662610100786, cgnt(0, 1, 1), 1e-12);
  EXPECT_NEAR(0.28209479177387842, cgnt(1, 1, 1), 1e-12);
  EXPECT_NEAR(0.0000000000000000, cgnt(2, 2, 2), 1e-12);
}