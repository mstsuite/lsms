//
// Created by F.Moitzi on 18.09.2021.
//

#undef NDEBUG

#include <assert.h>

#include <cstring>
#include <cstdlib>
#include <vector>
#include <cmath>

#include "integrator.hpp"
#include "poisson.hpp"

#include "Misc/integrateOneDim.hpp"


using namespace lsms;

template<typename T>
bool approx_equal(T x, T y, T epsilon) {
  return fabs(x - y) / max(fabs(x), fabs(y)) <= epsilon;
}

template<typename T>
T relative_diff(T ref, T val) {
  return std::fabs(ref - val) / std::fabs(ref);
}

int main(int argc, char *argv[]) {

  constexpr auto number_of_points = 800;

  std::vector<double> radial_mesh(number_of_points, 0.0);
  std::vector<double> radial_mesh_deriv(number_of_points, 0.0);
  std::vector<double> function(number_of_points, 0.0);
  std::vector<double> density(number_of_points, 0.0);
  std::vector<double> vhartreederiv(number_of_points, 0.0);
  std::vector<double> vhartree(number_of_points, 0.0);


  auto r0 = 0.00001;
  auto rmax = 2.0;

  auto h = std::log(rmax / r0) / number_of_points;


  for (auto i = 0; i < number_of_points; i++) {
    radial_mesh[i] = r0 * exp(i * h);
    radial_mesh_deriv[i] = r0 * exp(i * h) * h;
  }

  /*
   * 1. Test
   */

  for (auto i = 0; i < number_of_points; i++) {
    density[i] = 4 * M_PI * radial_mesh[i] * radial_mesh[i] * std::exp(-0.01 * radial_mesh[i]);
  }


  radial_poisson(vhartree, vhartreederiv, radial_mesh, radial_mesh_deriv, density, number_of_points);


  return EXIT_SUCCESS;
}