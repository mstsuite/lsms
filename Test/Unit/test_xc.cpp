//
// Created by F.Moitzi on 17.04.2022.
//

/**
 * Test suite for the calculation of the XC functional
 *
 * Compare XC potential and energy
 *
 */

#include <gtest/gtest.h>
#include <xc.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>

#include "Matrix.hpp"
#include "XCBase.hpp"
#include "XCLibxc.hpp"
#include "helpers.hpp"
#include "xc.hpp"

/**
 * Interface to built-in XC functionals
 */
extern "C" {

double alpha2_(double *rs, double *dz, double *sp, int *iexch, double *exchg);
}

/**
 *
 * The even more important problem is why does this not work for
 *
 */

namespace xc_tests {

double rsFromRho(double rho) { return std::cbrt(3.0 / (4.0 * M_PI * rho)); }

class XCLDATest : public ::testing::Test {
 protected:
  const int N = 10;
  const int nSpin = 1;
  std::vector<double> rho;
  Matrix<double> rhoIn;
  std::vector<double> rmesh;
  std::vector<double> drmesh;

  std::vector<double> vxc_ref;
  std::vector<double> exc_ref;

  XCLDATest() {
    vxc_ref = {-0.263791131116203059, -2.288466541003923993,
               -3.550395541899215246, -4.604164434312637333,
               -5.542561711608996688, -6.403761238593741467,
               -7.208160038861961105, -7.968265974292418896,
               -8.692444477225604871, -9.386637207703939723};

    exc_ref = {-0.203140250649610887, -1.740464242092732761,
               -2.692226997951151723, -3.485854482048946767,
               -4.192069046232443696, -4.839882036384272723,
               -5.444767923445450464, -6.016204512267600890,
               -6.560524936169307431, -7.082224164606909511};
  }

  void SetUp() override {
    rhoIn = Matrix<double>(N, nSpin);
    rho = std::vector<double>(N);
    rmesh = lsms_helper::linspace<double>(0.01, 10, N);
    drmesh = rmesh;
    for (int i = 0; i < N; i++) {
      drmesh[i] = rmesh[1] - rmesh[0];
      rho[i] = rmesh[i] * rmesh[i] + 0.001;
      rhoIn(i, 0) = rho[i] * rmesh[i] * rmesh[i] * 4 * M_PI;
    }
  }
};

TEST_F(XCLDATest, XCInterface) {
  Matrix<double> xcEnergyOut(N, 1);
  Matrix<double> xcPotOut(N, 1);

  lsms::XCLibxc xc(1, {1, 1, 7});

  // Check functionals

  auto &functionals = xc.get_functionals();

  ASSERT_STREQ(functionals[0].get_functional().info->name, "Slater exchange");
  ASSERT_STREQ(functionals[1].get_functional().info->name,
               "Vosko, Wilk & Nusair (VWN5)");

  xc.evaluate(rmesh, drmesh, rhoIn, N, xcEnergyOut, xcPotOut);

  for (int i = 0; i < N; i++) {
    EXPECT_NEAR(xcEnergyOut(i, 0), exc_ref[i], 1e-7);
    EXPECT_NEAR(xcPotOut(i, 0), vxc_ref[i], 1e-7);
  }
}

TEST_F(XCLDATest, BuiltInVWN) {
  std::vector<double> vxc(N);
  std::vector<double> exc(N);

  double dz = 0.0;
  double sp = 0.0;
  int iexch = 2;

  for (int i = 0; i < N; i++) {
    auto rs = rsFromRho(rho[i]);
    vxc[i] = alpha2_(&rs, &dz, &sp, &iexch, &exc[i]);
  }

  for (int i = 0; i < N; i++) {
    EXPECT_NEAR(exc[i], exc_ref[i], 1e-7);
    EXPECT_NEAR(vxc[i], vxc_ref[i], 1e-7);
  }
}

TEST_F(XCLDATest, LibxcVWN) {
  std::vector<double> vxc(N);
  std::vector<double> exc(N);

  const int libxc_lda_x = 1;
  const int libxc_lda_BH = 7;

  std::vector<double> vxc_libxc_x(N, 0.0);
  std::vector<double> exc_libxc_x(N, 0.0);

  std::vector<double> vxc_libxc_c(N, 0.0);
  std::vector<double> exc_libxc_c(N, 0.0);

  xc_func_type xFunctional, cFunctional;

  EXPECT_EQ(xc_func_init(&xFunctional, libxc_lda_x, XC_UNPOLARIZED), 0);
  EXPECT_EQ(xc_func_init(&cFunctional, libxc_lda_BH, XC_UNPOLARIZED), 0);

  xc_lda_exc_vxc(&xFunctional, N, rho.data(), exc_libxc_x.data(),
                 vxc_libxc_x.data());

  xc_lda_exc_vxc(&cFunctional, N, rho.data(), exc_libxc_c.data(),
                 vxc_libxc_c.data());

  std::transform(exc_libxc_x.begin(), exc_libxc_x.end(), exc_libxc_c.begin(),
                 exc.begin(), std::plus<>());

  std::transform(vxc_libxc_x.begin(), vxc_libxc_x.end(), vxc_libxc_c.begin(),
                 vxc.begin(), std::plus<>());

  std::transform(vxc.begin(), vxc.end(), vxc.begin(), [](auto number) {
    return 2.0 * number;
  });  // negates elements correctly

  std::transform(exc.begin(), exc.end(), exc.begin(), [](auto number) {
    return 2.0 * number;
  });  // negates elements correctly

  for (int i = 0; i < N; i++) {
    EXPECT_NEAR(exc[i], exc_ref[i], 1e-7);
    EXPECT_NEAR(vxc[i], vxc_ref[i], 1e-7);
  }
}

TEST_F(XCLDATest, ReferVWN) {
  std::vector<double> vxc(N);
  std::vector<double> exc(N);

  bool relat = false;
  double c_light = 137.0359895;

  for (int i = 0; i < N; i++) {
    getvxc_scalar(&rho[i], &relat, &c_light, &exc[i], &vxc[i]);
  }

  std::transform(exc.begin(), exc.end(), exc.begin(), [](auto number) {
    return 2.0 * number;
  });  // negates elements correctly

  std::transform(vxc.begin(), vxc.end(), vxc.begin(), [](auto number) {
    return 2.0 * number;
  });  // negates elements correctly

  for (int i = 0; i < N; i++) {
    EXPECT_NEAR(exc[i], exc_ref[i], 1e-7);
    EXPECT_NEAR(vxc[i], vxc_ref[i], 1e-7);
  }
}

class XCGGATest : public ::testing::Test {
 protected:
  const int nSpin = 1;

  Matrix<double> rhoIn;

  std::vector<double> rho;
  std::vector<double> rmesh;
  std::vector<double> drmesh;

  std::vector<double> vxc_ref;
  std::vector<double> exc_ref;

  XCGGATest() {
    std::string filename("result_pbe_libxc.out");
    std::ifstream data(filename.c_str());

    // Read in the test file
    if (data.is_open()) {
      double r;
      double rrho;
      double v;
      double e;
      std::vector<double> line;
      while (data >> v >> e >> r >> rrho) {
        rmesh.push_back(r);
        rho.push_back(rrho);
        vxc_ref.push_back(v);
        exc_ref.push_back(e);
      }
    }

    rhoIn = Matrix<double>(rho.size(), nSpin);

    drmesh = rmesh;

    for (int ir = 0; ir < rho.size(); ir++) {
      rhoIn(ir, 0) = rho[ir];
      drmesh[ir] = rmesh[ir] * 0.005;
    }
  }
};

TEST_F(XCGGATest, XCInterface) {
  Matrix<double> xcEnergyOut(rho.size(), 1);
  Matrix<double> xcPotOut(rho.size(), 1);

  lsms::XCLibxc xc(1, {1, 101, 130});

  // Check functionals
  auto &functionals = xc.get_functionals();

  ASSERT_STREQ(functionals[0].get_functional().info->name,
               "Perdew, Burke & Ernzerhof");
  ASSERT_STREQ(functionals[1].get_functional().info->name,
               "Perdew, Burke & Ernzerhof");

  int N = 1000;

  xc.evaluate(rmesh, drmesh, rhoIn, N, xcEnergyOut, xcPotOut);

  /*
   * Check if energy density is correctly calculated
   */
  for (int i = 0; i < N; i++) {
    ASSERT_NEAR(xcEnergyOut(i, 0), exc_ref[i], 5e-7);
    // std::printf("%30.20f, %30.20f\n", xcEnergyOut(i, 0), exc_ref[i]);
  }

  for (int i = 0; i < N; i++) {
    ASSERT_NEAR(xcPotOut(i, 0), vxc_ref[i], 2e-3);
    // std::printf("%30.20f, %30.20f\n", xcPotOut(i, 0), vxc_ref[i]);
  }
}

}  // namespace xc_tests