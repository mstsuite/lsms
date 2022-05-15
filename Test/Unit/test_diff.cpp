//
// Created by F.Moitzi on 15.05.2022.
//

#include <gtest/gtest.h>
#include <cmath>

#include "Misc/diff.hpp"

namespace diff_tests {


  TEST(DiffTests, SimpleMesh) {

    constexpr auto N = 1000;
    std::vector<double> fun(N);
    std::vector<double> ref(N);
    std::vector<double> x(N);

    for (auto ir{0}; ir < N; ir++) {
      x[ir] = ir;
      fun[ir] = x[ir] * x[ir] + 1.0;
      ref[ir] = 2.0 * x[ir];
    }

    auto res = lsms::derivative<double>(fun, N);

    for (auto ir{0}; ir < N; ir++) {
      EXPECT_NEAR(ref[ir], res[ir], 10e-6);
    }

  }

  TEST(DiffTests, ExponetialMesh) {

    constexpr auto N = 1000;
    std::vector<double> fun(N);
    std::vector<double> ref(N);
    std::vector<double> r(N);

    for (auto ir{0}; ir < N; ir++) {
      r[ir] = 0.01 * std::exp(ir * 0.01);
      fun[ir] = r[ir] * r[ir] + 1.0;
      ref[ir] = 2.0 * r[ir];
    }

    auto res = lsms::derivative<double>(fun, N);

    for (auto ir{0}; ir < N; ir++) {
      res[ir] = res[ir] / (r[ir] * 0.01);
    }

    for (auto ir{0}; ir < N; ir++) {
      EXPECT_NEAR(ref[ir], res[ir], 10e-5);
    }

  }

}