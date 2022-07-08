//
// Created by F.Moitzi on 15.05.2022.
//
#include <gtest/gtest.h>

#include <fstream>
#include <iostream>

#include "Misc/integrateOneDim.hpp"
#include "Misc/integrator.hpp"

/**
 * Test suite for the integration
 *
 * This tests are for the integration routines mostly used in
 * `localTotalEnergy.cpp`. They also act as a benchmarks.
 *
 * The input densities correspond to 3 integral that were obtain from
 * 3 consecutive runs with converged potentials. It shows that the
 * integration routines of LSMS have a large numerical error.
 *
 */
namespace integrate_energy_tests {

struct IntegrateEnergyTestFixture : public testing::TestWithParam<std::string> {
};

TEST_P(IntegrateEnergyTestFixture, Ezrho) {
  std::string filename = GetParam();

  std::vector<double> rmesh;
  std::vector<double> integrand;

  std::ifstream data(filename.c_str());

  if (data.is_open()) {
    double r;
    double val;
    std::vector<double> line;
    while (data >> val >> r) {
      rmesh.push_back(r);
      integrand.push_back(val);
    }
  }

  std::vector<double> integral = integrand;

  int ir_max = 500;
  double rSphere = rmesh[ir_max - 100];

  // Index of element that is less then `r_sphere`
  auto iter = std::find_if(rmesh.rbegin(), rmesh.rend(), [rSphere](double d) {
    if (std::fabs(rSphere - d) < 1e-5) {
      return true;
    } else {
      return false;
    }
  });

  auto end_point = std::distance(rmesh.begin(), iter.base());

  auto result_11 = integrateOneDim(rmesh, integrand, integral, rSphere);
  auto result_12 = lsms::radialIntegral(integrand, rmesh, end_point);

  std::printf("%20.10f\n", result_11);
  std::printf("%20.10f\n", result_12);

  EXPECT_NEAR(result_11, 746.9981173558, 10e-2);
  EXPECT_NEAR(result_12, 746.9981173558, 10e-8);
}

INSTANTIATE_TEST_SUITE_P(IntegrateEnergyTests, IntegrateEnergyTestFixture,
                        ::testing::Values("ezrho1.out", "ezrho2.out",
                                          "ezrho3.out"));

}  // namespace integrate_energy_tests
