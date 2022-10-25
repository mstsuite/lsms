//
// Created by F.Moitzi on 18.09.2021.
//

#include <gtest/gtest.h>

#undef NDEBUG

#include <cmath>
#include <vector>

#include "Misc/integrateOneDim.hpp"
#include "Misc/integrator.hpp"
#include "Misc/poisson.hpp"
#include "accel_common.hpp"

namespace poisson_tests {

  using namespace lsms;

  template<typename T>
  T relative_diff(T ref, T val) {
    return std::fabs(ref - val) / std::fabs(ref);
  }

/**
 * Tests for the radial poission equation
 *
 * This equation is a Euler-Cauchy equation that can be solved analytically for
 * certain cases
 *
 */
  TEST(IntegrationPoisson, CompareDirect) {
    constexpr auto number_of_points = 800;

    std::vector<double> radial_mesh(number_of_points, 0.0);
    std::vector<double> radial_mesh_deriv(number_of_points, 0.0);
    std::vector<double> function(number_of_points, 0.0);
    std::vector<double> density(number_of_points, 0.0);

    std::vector<double> integral(number_of_points, 0.0);
    std::vector<double> integrand(number_of_points, 0.0);

    std::vector<double> reference(number_of_points, 0.0);

    std::vector<double> vhartreederiv(number_of_points, 0.0);
    std::vector<double> vhartree(number_of_points, 0.0);

    auto r0 = 0.00001;
    auto rmax = 2.0;

    auto h = std::log(rmax / r0) / (number_of_points - 1);

    for (auto i = 0; i < number_of_points; i++) {
      radial_mesh[i] = r0 * exp(i * h);
      radial_mesh_deriv[i] = r0 * exp(i * h) * h;
    }

    /*
     * 1. Test
     *
     * \rho(r) = 4 * pi * r^2
     *
     * Analytic reference solution:
     *
     * V(r) = - 2/3 * pi * ( x*x - 12)
     *
     */

    for (auto i = 0; i < number_of_points; i++) {
      density[i] = 4 * M_PI * radial_mesh[i] * radial_mesh[i];
      auto &r = radial_mesh[i];
      reference[i] = -M_PI * 2.0 / 3.0 * (r * r - 12);
    }

    radial_poisson(vhartree, vhartreederiv, radial_mesh, radial_mesh_deriv,
                   density, number_of_points);

    for (auto i = 0; i < number_of_points; i++) {
      ASSERT_TRUE(relative_diff(reference[i], vhartree[i]) < 2e-11);
    }

    /*
     * 2. Test
     *
     * \rho(r) = 4 * pi * r^2 * exp(-r)
     *
     * Analytic reference solution:
     *
     * V(r) = 4 * pi * (- e^(-x) * (x + 2) / x + 2/ x - 3/ e^(2))
     *
     */

    for (auto i = 0; i < number_of_points; i++) {
      auto &r = radial_mesh[i];
      density[i] = 4 * M_PI * r * r * exp(-r);
      reference[i] = 4 * M_PI * (-exp(-r) * (r + 2) / r + 2 / r - 3 / exp(2));
    }

    radial_poisson(vhartree, vhartreederiv, radial_mesh, radial_mesh_deriv,
                   density, number_of_points);

    for (auto i = 0; i < number_of_points; i++) {
      ASSERT_TRUE(relative_diff(reference[i], vhartree[i]) < 1e-7);
    }

    /*
     * 3. Test
     *
     * Test calculation of hartree energy
     *
     * \rho(r) = 4 * pi * r^2 * exp(-r)
     *
     * V(r) = 4 * pi * (- e^(-x) * (x + 2) / x + 2/ x - 3/ e^(2))
     *
     * e^−4 * (20 * e^4 − 192 * e^2 + 572 ) * pi^2
     *
     */

    for (auto i = 0; i < number_of_points; i++) {
      auto &r = radial_mesh[i];
      density[i] = 4 * M_PI * r * r * exp(-r);
      reference[i] = 4 * M_PI * (-exp(-r) * (r + 2) / r + 2 / r - 3 / exp(2));
    }

    for (auto i = 0; i < number_of_points; i++) {
      integrand[i] = density[i];
    }
    integrateOneDim(radial_mesh, integrand, integral,
                    radial_mesh[number_of_points - 1]);
    for (auto i = 0; i < number_of_points; i++) {
      integral[i] = 2.0 * integral[i] * density[i] / radial_mesh[i];
    }
    auto reference_energy = integrateOneDim(radial_mesh, integral, integrand,
                                            radial_mesh[number_of_points - 1]);

    radial_poisson(vhartree, vhartreederiv, radial_mesh, radial_mesh_deriv,
                   density, number_of_points);
    for (auto i = 0; i < number_of_points; i++) {
      integral[i] = vhartree[i] * density[i];
    }
    auto energy = radialIntegral(integral, radial_mesh, number_of_points);

    auto analytical_energy =
        exp(-4) * (20 * exp(4) - 192 * exp(2) + 572) * M_PI * M_PI;

    ASSERT_TRUE(relative_diff(analytical_energy, reference_energy) < 3e-6);
    ASSERT_TRUE(relative_diff(analytical_energy, energy) < 6e-8)
                  << "All tests have finished";
  }

  TEST(IntegrationPoisson, CompareExponetialMesh) {
    constexpr auto number_of_points = 800;

    std::vector<double> radial_mesh(number_of_points, 0.0);
    std::vector<double> radial_mesh_deriv(number_of_points, 0.0);
    std::vector<double> function(number_of_points, 0.0);
    std::vector<double> density(number_of_points, 0.0);

    std::vector<double> integral(number_of_points, 0.0);
    std::vector<double> integrand(number_of_points, 0.0);

    std::vector<double> reference(number_of_points, 0.0);

    std::vector<double> vhartreederiv(number_of_points, 0.0);
    std::vector<double> vhartree(number_of_points, 0.0);

    auto r0 = 0.00001;
    auto rmax = 2.0;

    auto h = std::log(rmax / r0) / (number_of_points - 1);

    for (auto i = 0; i < number_of_points; i++) {
      radial_mesh[i] = r0 * exp(i * h);
    }

    /*
     * 1. Test
     *
     * \rho(r) = 4 * pi * r^2
     *
     * Analytic reference solution:
     *
     * V(r) = - 2/3 * pi * ( x*x - 12)
     *
     */

    for (auto i = 0; i < number_of_points; i++) {
      density[i] = 4 * M_PI * radial_mesh[i] * radial_mesh[i];
      auto &r = radial_mesh[i];
      reference[i] = -M_PI * 2.0 / 3.0 * (r * r - 12);
    }

    radial_poisson(vhartree, vhartreederiv, radial_mesh, h, density,
                   number_of_points);

    for (auto i = 0; i < number_of_points; i++) {
      ASSERT_TRUE(relative_diff(reference[i], vhartree[i]) < 2e-11);
    }

    /*
     * 2. Test
     *
     * \rho(r) = 4 * pi * r^2 * exp(-r)
     *
     * Analytic reference solution:
     *
     * V(r) = 4 * pi * (- e^(-x) * (x + 2) / x + 2/ x - 3/ e^(2))
     *
     */

    for (auto i = 0; i < number_of_points; i++) {
      auto &r = radial_mesh[i];
      density[i] = 4 * M_PI * r * r * exp(-r);
      reference[i] = 4 * M_PI * (-exp(-r) * (r + 2) / r + 2 / r - 3 / exp(2));
    }

    radial_poisson(vhartree, vhartreederiv, radial_mesh, h, density,
                   number_of_points);

    for (auto i = 0; i < number_of_points; i++) {
      ASSERT_TRUE(relative_diff(reference[i], vhartree[i]) < 1e-7);
    }
  }

  TEST(IntegrationPoisson, CompareExponetialMeshCorrection) {
    constexpr auto number_of_points = 400;

    std::vector<double> radial_mesh(number_of_points, 0.0);
    std::vector<double> radial_mesh_deriv(number_of_points, 0.0);
    std::vector<double> function(number_of_points, 0.0);
    std::vector<double> density(number_of_points, 0.0);

    std::vector<double> integral(number_of_points, 0.0);
    std::vector<double> integrand(number_of_points, 0.0);

    std::vector<double> reference(number_of_points, 0.0);

    std::vector<double> vhartreederiv(number_of_points, 0.0);
    std::vector<double> vhartree(number_of_points, 0.0);

    auto r0 = 0.001;
    auto rmax = 2.0;

    auto h = std::log(rmax / r0) / (number_of_points - 1);

    for (auto i = 0; i < number_of_points; i++) {
      radial_mesh[i] = r0 * exp(i * h);
    }

    /*
     * 1. Test
     *
     * \rho(r) = 4 * pi * r^2
     *
     * Analytic reference solution:
     *
     * V(r) = - 2/3 * pi * ( x*x - 12)
     *
     */

    for (auto i = 0; i < number_of_points; i++) {
      density[i] = 4 * M_PI * radial_mesh[i] * radial_mesh[i];
      auto &r = radial_mesh[i];
      reference[i] = -M_PI * 2.0 / 3.0 * (r * r - 12);
    }

    radial_poisson(vhartree, vhartreederiv, radial_mesh, h, density,
                   number_of_points);

    for (auto i = 0; i < number_of_points; i++) {
      std::cout << reference[i] << " " << vhartree[i] << relative_diff(reference[i], vhartree[i]) << std::endl;

      //ASSERT_TRUE(relative_diff(reference[i], vhartree[i]) < 2e-11);
    }


  }


}  // namespace poisson_tests