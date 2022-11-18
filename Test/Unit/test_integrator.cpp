//
// Created by F.Moitzi on 18.09.2021.
//

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "Misc/integrateOneDim.hpp"
#include "accel_common.hpp"
#include "integrator.hpp"

namespace integrator_tests {

using namespace lsms;

template <typename T>
T relative_diff(T ref, T val) {
  return std::fabs(ref - val) / std::fabs(ref);
}

TEST(IntegrationTests, CompareRoutines) {
  constexpr auto number_of_points = 800;

  std::vector<double> radial_mesh(number_of_points, 0.0);
  std::vector<double> radial_mesh_deriv(number_of_points, 0.0);
  std::vector<double> function(number_of_points, 0.0);

  std::vector<double> integrand(number_of_points, 0.0);
  std::vector<double> integral(number_of_points, 0.0);

  auto r0 = 0.00001;
  auto h = 0.0155;

  for (auto i = 0; i < number_of_points; i++) {
    radial_mesh[i] = r0 * exp(i * h);
    radial_mesh_deriv[i] = r0 * exp(i * h) * h;
  }

  auto rSphere = radial_mesh[number_of_points - 1];

  /*
   * 1. Function test
   */

  auto func = [](auto x) { return (20 * x * x * exp(-2.0 * x)); };

  for (auto i = 0; i < number_of_points; i++) {
    function[i] = func(radial_mesh[i]);
  }

  constexpr auto reference_with_zero = 4.278412171945034;
  constexpr auto reference_without_zero = 4.278412171945027;

  auto result_01 = integrateOneDim(radial_mesh, function, integral, rSphere);
  auto result_02 = radialIntegral(function, radial_mesh, number_of_points);
  auto result_03 = radialIntegral(function, radial_mesh, number_of_points);
  auto result_04 =
      radialIntegralDerivMesh(function, radial_mesh_deriv, number_of_points);

#ifdef PRINT
  std::printf("%30.15f\n", relative_diff(reference_without_zero, result_01));
  std::printf("%30.15f\n", relative_diff(reference_without_zero, result_02));
  std::printf("%30.15f\n", relative_diff(reference_without_zero, result_03));
  std::printf("%30.15f\n", relative_diff(reference_without_zero, result_04));
#endif

  ASSERT_TRUE(relative_diff(reference_without_zero, result_01) < 1e-6);
  ASSERT_TRUE(relative_diff(reference_without_zero, result_02) < 1e-8);
  ASSERT_TRUE(relative_diff(reference_without_zero, result_03) < 1e-8);
  ASSERT_TRUE(relative_diff(reference_without_zero, result_04) < 1e-6);

  /*
   * 2. Function test
   */

  auto func_2 = [](auto x) { return (0.02 * x * x * x - x * x + 20 - 40 * x); };

  constexpr auto reference_2_without_zero = -70.88417469817453;

  for (auto i = 0; i < number_of_points; i++) {
    function[i] = func_2(radial_mesh[i]);
  }

  auto result_11 = integrateOneDim(radial_mesh, function, integral, rSphere);
  auto result_12 = radialIntegral(function, radial_mesh, number_of_points);
  auto result_13 = radialIntegral(function, radial_mesh, number_of_points);
  auto result_14 =
      radialIntegralDerivMesh(function, radial_mesh_deriv, number_of_points);

#ifdef PRINT
  std::printf("%30.15f\n", relative_diff(reference_2_without_zero, result_11));
  std::printf("%30.15f\n", relative_diff(reference_2_without_zero, result_12));
  std::printf("%30.15f\n", relative_diff(reference_2_without_zero, result_13));
  std::printf("%30.15f\n", relative_diff(reference_2_without_zero, result_14));
#endif

  ASSERT_TRUE(relative_diff(reference_2_without_zero, result_11) < 1e-6);
  ASSERT_TRUE(relative_diff(reference_2_without_zero, result_12) < 1e-9);
  ASSERT_TRUE(relative_diff(reference_2_without_zero, result_13) < 1e-9);
  ASSERT_TRUE(relative_diff(reference_2_without_zero, result_14) < 1e-5);

  /*
   * 3. Function test
   */

  constexpr auto reference_3_without_zero = 2.691512828088295;

  auto func_3 = [](auto x) { return (sin(x) - x * x * cos(x)); };

  for (auto i = 0; i < number_of_points; i++) {
    function[i] = func_3(radial_mesh[i]);
  }

  auto result_21 = integrateOneDim(radial_mesh, function, integral, rSphere);
  auto result_22 = radialIntegral(function, radial_mesh, number_of_points);
  auto result_23 = radialIntegral(function, radial_mesh, number_of_points);
  auto result_24 =
      radialIntegralDerivMesh(function, radial_mesh_deriv, number_of_points);

#ifdef PRINT
  std::printf("%30.15f\n", relative_diff(reference_3_without_zero, result_21));
  std::printf("%30.15f\n", relative_diff(reference_3_without_zero, result_22));
  std::printf("%30.15f\n", relative_diff(reference_3_without_zero, result_23));
  std::printf("%30.15f\n", relative_diff(reference_3_without_zero, result_24));
#endif

  ASSERT_TRUE(relative_diff(reference_3_without_zero, result_21) < 1e-5);
  ASSERT_TRUE(relative_diff(reference_3_without_zero, result_22) < 1e-6);
  ASSERT_TRUE(relative_diff(reference_3_without_zero, result_23) < 1e-6);
  ASSERT_TRUE(relative_diff(reference_3_without_zero, result_24) < 1e-4);
}

}  // namespace integrator_tests