//
// Created by F.Moitzi on 18.09.2021.
//

#include <gtest/gtest.h>

#include "GauntFactor.hpp"
#include "Indices.hpp"
#include "Coeficients.hpp"

TEST(MathTests, GauntCoeffiencents1) {

  int lmax = 20;

  lsms::math::GauntFactor gaunt_factor(lmax);

  // First tests
  {
    auto l1 = 2;
    auto l2 = 2;
    auto l3 = 0;

    auto m1 = 0;
    auto m2 = 0;
    auto m3 = 0;

    auto k1 = l1 * l1 + l1 - m1;
    auto k2 = l2 * l2 + l2 - m2;
    auto k3 = l3 * l3 + l3 - m3;

    EXPECT_NEAR(0.2820947917738781434740397257803862929220253146644994284220428609,
                gaunt_factor.table(k1, k2, k3), 1e-8);
  }

  // Second tests
  {
    auto l1 = 4;
    auto l2 = 2;
    auto l3 = 4;

    auto m1 = 1;
    auto m2 = 2;
    auto m3 = 1;

    auto k1 = l1 * l1 + l1 - m1;
    auto k2 = l2 * l2 + l2 - m2;
    auto k3 = l3 * l3 + l3 - m3;

    EXPECT_NEAR(0.2006619231289296521258144079540202145901815398377730367967227404,
                gaunt_factor.table(k1, k2, k3), 1e-8);
  }

  // Third tests
  {

    auto l1 = 6;
    auto l2 = 4;
    auto l3 = 4;

    auto m1 = 2;
    auto m2 = 3;
    auto m3 = 1;

    auto k1 = l1 * l1 + l1 - m1;
    auto k2 = l2 * l2 + l2 - m2;
    auto k3 = l3 * l3 + l3 - m3;

    EXPECT_NEAR(0.1652827715004524673443686414985987229164078374402080257136786467,
                gaunt_factor.table(k1, k2, k3), 1e-8);
  }

}

TEST(MathTests, GauntCoeffiencents2) {

  int lmax = 3;

  lsms::math::GauntFactor gaunt_factor(lmax);

  {
    auto l1 = 3;
    auto l2 = 3;
    auto l3 = 3;

    auto m1 = 0;
    auto m2 = 0;
    auto m3 = 0;

    auto k1 = l1 * l1 + l1 - m1;
    auto k2 = l2 * l2 + l2 - m2;
    auto k3 = l3 * l3 + l3 - m3;

    gaunt_factor.table(k1, k2, k3);
  }


}

TEST(MathTests, GauntCoeffiencents3) {

  int lmax = 6;

  SphericalHarmonicsCoeficients s;
  s.init(lmax);

  AngularMomentumIndices a;
  a.init(lmax);

  std::vector<Real> tg, wg;
  Matrix<Real> plmg;
  Array3d<Real> cgnt;

  tg.resize(2 * (2 * lmax + 1));
  wg.resize(2 * (2 * lmax + 1));
  cgnt.resize(lmax + 1, (lmax + 1) * (lmax + 1), (lmax + 1) * (lmax + 1));

  plmg.resize(((2 * lmax + 1) * (2 * lmax + 2)) / 2, 2 * lmax + 1);

  int iprint = -1;
  char *istop = 0;

  cgaunt_(&lmax, &s.clm[0], &plmg(0, 0), &tg[0], &wg[0],
          &cgnt(0, 0, 0),
          &a.lofk[0], &a.mofk[0],
          &iprint, istop);


  // First tests
  {
    auto l1 = 0;
    auto l2 = 0;
    auto l3 = 0;

    auto m1 = 0;
    auto m2 = 0;

    // Not necessary because it is zero otherwise
    // auto m3 = 0;

    auto k1 = l1 * l1 + l1 - m1;
    auto k2 = l2 * l2 + l2 - m2;

    auto w3 = l3 / 2;

    EXPECT_NEAR(0.2820947917738781434740397257803862929220253146644994284220428609,
                cgnt(w3, k1, k2), 1e-8);
  }

  // Second tests
  {
    auto l1 = 4;
    auto l2 = 2;
    auto l3 = 4;

    auto m1 = 1;
    auto m2 = 2;

    // Not necessary because it is zero otherwise
    // auto m3 = 0;

    auto k1 = l1 * l1 + l1 - m1;
    auto k2 = l2 * l2 + l2 - m2;
    auto w3 = l3 / 2;

    EXPECT_NEAR(0.2006619231289296521258144079540202145901815398377730367967227404,
                cgnt(w3, k1, k2), 1e-8);
  }


  // Third tests
  {

    auto l1 = 6;
    auto l2 = 4;
    auto l3 = 4;

    auto m1 = 2;
    auto m2 = 3;

    auto k1 = l1 * l1 + l1 - m1;
    auto k2 = l2 * l2 + l2 - m2;
    auto w3 = l3 / 2;

    EXPECT_NEAR(0.1652827715004524673443686414985987229164078374402080257136786467,
                cgnt(w3, k1, k2), 1e-8);
  }


}