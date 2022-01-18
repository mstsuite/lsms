//
// Created by F.Moitzi on 23.09.2021.
//

#ifndef LSMS_INTEGRATOR_HPP
#define LSMS_INTEGRATOR_HPP

#include <algorithm>
#include <functional>

#include "utils.hpp"

namespace lsms {

  template<typename T, typename U, std::enable_if_t<std::is_integral<U>::value, bool> = true>
  std::vector<T> diff(const std::vector<T> &x, U length) {

    std::vector<T> result(length);

    for (U i = 0; i < length - 1; i++) {
      result[i] = x[i + 1] - x[i];
    }

    return result;
  }


  template<class T1, class T2>
  T2 simpson_nonuniform(std::vector<T1> x, std::vector<T2> f, std::size_t length) {

    auto N = length - 1;
    auto h = diff(x, length);

    T2 result {0.0};

    for (auto i = 1; i < N; i += 2) {

      auto hph = h[i] + h[i - 1];

      result += f[i] * (lsms::math::pow(h[i], 3) + lsms::math::pow(h[i - 1], 3)
                        + 3.0 * h[i] * h[i - 1] * hph)
                / (6 * h[i] * h[i - 1]);

      result += f[i - 1] * (2.0 * lsms::math::pow(h[i - 1], 3) - lsms::math::pow(h[i], 3)
                            + 3.0 * h[i] * lsms::math::pow(h[i - 1], 2))
                / (6 * h[i - 1] * hph);

      result += f[i + 1] * (2.0 * lsms::math::pow(h[i], 3) - lsms::math::pow(h[i - 1], 3)
                            + 3.0 * h[i - 1] * lsms::math::pow(h[i], 2))
                / (6 * h[i] * hph);

    }


    if ((N + 1) % 2 == 0) {

      if (N > 1) {

        result += f[N] * (2 * lsms::math::pow(h[N - 1], 2)
                          + 3.0 * h[N - 2] * h[N - 1])
                  / (6 * (h[N - 2] + h[N - 1]));

        result += f[N - 1] * (lsms::math::pow(h[N - 1], 2)
                              + 3.0 * h[N - 1] * h[N - 2])
                  / (6 * h[N - 2]);

        result -= f[N - 2] * lsms::math::pow(h[N - 1], 3)
                  / (6 * h[N - 2] * (h[N - 2] + h[N - 1]));

      } else {

        result += 0.5 * (f[N] + f[N - 1]) * h[N - 1];

      }
    }


    return result;
  }

  /**
   * Radial integrator for any mesh
   *
   * This method doesn't use any interpolations for the last point and assumes an alight mesh.
   * The mesh should be align to according to the `generateRadialMesh` method
   */
  template<typename T>
  T radialIntegral(const std::vector<T> &integrand,
                   const std::vector<T> &radial_mesh,
                   const T r_sphere) {


    // Index of element that is less then `r_sphere`
    auto iter = std::upper_bound(radial_mesh.rbegin(), radial_mesh.rend(), r_sphere, std::greater<T>());

    // Index of element that has radius greater then equal `r_sphere`
    auto i_after_rs = std::distance(radial_mesh.begin(), iter.base());

    auto end_point = i_after_rs + 1;

    return simpson_nonuniform(radial_mesh, integrand, end_point);
  }

}


#endif //LSMS_INTEGRATOR_HPP
