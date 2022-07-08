//
// Created by F.Moitzi on 23.09.2021.
//

#ifndef LSMS_INTEGRATOR_HPP
#define LSMS_INTEGRATOR_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <type_traits>
#include <vector>

#include "utils.hpp"

namespace lsms {

template <typename T, typename U,
          std::enable_if_t<std::is_integral<U>::value, bool> = true>
std::vector<T> diff(const std::vector<T> &x, U length) {
  std::vector<T> result(length);

  for (U i = 0; i < length - 1; i++) {
    result[i] = x[i + 1] - x[i];
  }

  return result;
}

template <typename T, typename U,
          std::enable_if_t<std::is_integral<U>::value, bool> = true>
T simpson_nonuniform(std::vector<T> x, std::vector<T> f, U length) {
  auto N = length - 1;
  auto h = diff(x, length);

  auto result = 0.0;

  for (U i = 1; i < N; i += 2) {
    auto hph = h[i] + h[i - 1];

    result += f[i] *
              (lsms::pow(h[i], 3) + lsms::pow(h[i - 1], 3) +
               3.0 * h[i] * h[i - 1] * hph) /
              (6 * h[i] * h[i - 1]);

    result += f[i - 1] *
              (2.0 * lsms::pow(h[i - 1], 3) - lsms::pow(h[i], 3) +
               3.0 * h[i] * lsms::pow(h[i - 1], 2)) /
              (6 * h[i - 1] * hph);

    result += f[i + 1] *
              (2.0 * lsms::pow(h[i], 3) - lsms::pow(h[i - 1], 3) +
               3.0 * h[i - 1] * lsms::pow(h[i], 2)) /
              (6 * h[i] * hph);
  }

  if ((N + 1) % 2 == 0) {
    if (N > 1) {
      result += f[N] *
                (2 * lsms::pow(h[N - 1], 2) + 3.0 * h[N - 2] * h[N - 1]) /
                (6 * (h[N - 2] + h[N - 1]));

      result += f[N - 1] *
                (lsms::pow(h[N - 1], 2) + 3.0 * h[N - 1] * h[N - 2]) /
                (6 * h[N - 2]);

      result -= f[N - 2] * lsms::pow(h[N - 1], 3) /
                (6 * h[N - 2] * (h[N - 2] + h[N - 1]));

    } else {
      result += 0.5 * (f[N] + f[N - 1]) * h[N - 1];
    }
  }

  return result;
}

/**
 * Radial integrator for any mesh
 */
template <typename T>
T radialIntegral(const std::vector<T> &integrand,
                 const std::vector<T> &radial_mesh, std::size_t end) {
  auto result = simpson_nonuniform(radial_mesh, integrand, end);

  return result;
}

/**
 * Radial integrator for mesh with an analytic derivative
 *
 * \f[ r = r_0 * exp(i * h) ]\f
 *
 * dr / di = r * h
 *
 */
template <typename T>
T radialIntegralDerivMesh(const std::vector<T> &val,
                          const std::vector<T> &radial_mesh_deriv, std::size_t end) {
  T intval = 0.0;

  auto idx = end;
  for (idx = 0; idx < end - 2; idx += 2) {
    intval += 1.0 / 3.0 *
              (val[idx] * radial_mesh_deriv[idx] +
               4 * val[idx + 1] * radial_mesh_deriv[idx + 1] +
               val[idx + 2] * radial_mesh_deriv[idx + 2]);
  };

  if ((end % 2) == 0) {
    intval += 0.5 * (val[idx] * radial_mesh_deriv[idx] +
                     val[idx + 1] * radial_mesh_deriv[idx + 1]);
  };
  return intval;
}

}  // namespace lsms

#endif  // LSMS_INTEGRATOR_HPP
