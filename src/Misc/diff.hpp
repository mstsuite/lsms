//
// Created by F.Moitzi on 15.05.2022.
//

#ifndef LSMS_DIFF_HPP
#define LSMS_DIFF_HPP

#include <cstdlib>
#include <vector>

#include "Matrix.hpp"

namespace lsms {

template <typename T>
std::vector<T> derivative(const T *fun, std::size_t endpoint) {
  std::vector<T> dfun(endpoint);

  // 6-th order forward difference for the first two points
  dfun[0] = (-49.0 / 20.0 * fun[0] + 6.0 * fun[1] - 15.0 / 2.0 * fun[2] +
             20.0 / 3.0 * fun[3] - 15.0 / 4.0 * fun[4] + 6.0 / 5.0 * fun[5] -
             1.0 / 6.0 * fun[6]);

  dfun[1] = (-49.0 / 20.0 * fun[1] + 6.0 * fun[2] - 15.0 / 2.0 * fun[3] +
             20.0 / 3.0 * fun[4] - 15.0 / 4.0 * fun[5] + 6.0 / 5.0 * fun[6] -
             1.0 / 6.0 * fun[7]);

  // 4-th order central difference
  for (int ir = 2; ir < endpoint - 2; ir++) {
    dfun[ir] =
        (fun[ir - 2] - 8.0 * fun[ir - 1] + 8.0 * fun[ir + 1] - fun[ir + 2]) /
        12.0;
  }

  // 6-th order backward difference for the last two points
  dfun[endpoint - 2] =
      -(-49.0 / 20.0 * fun[endpoint - 2] + 6.0 * fun[endpoint - 3] -
        15.0 / 2.0 * fun[endpoint - 4] + 20.0 / 3.0 * fun[endpoint - 5] -
        15.0 / 4.0 * fun[endpoint - 6] + 6.0 / 5.0 * fun[endpoint - 7] -
        1.0 / 6.0 * fun[endpoint - 8]);

  dfun[endpoint - 1] =
      -(-49.0 / 20.0 * fun[endpoint - 1] + 6.0 * fun[endpoint - 2] -
        15.0 / 2.0 * fun[endpoint - 3] + 20.0 / 3.0 * fun[endpoint - 4] -
        15.0 / 4.0 * fun[endpoint - 5] + 6.0 / 5.0 * fun[endpoint - 6] -
        1.0 / 6.0 * fun[endpoint - 7]);

  return dfun;
}

template <typename T>
std::vector<T> derivative(const std::vector<T> &fun, std::size_t endpoint) {
  auto dfun = fun;

  // 6-th order forward difference for the first two points
  dfun[0] = (-49.0 / 20.0 * fun[0] + 6.0 * fun[1] - 15.0 / 2.0 * fun[2] +
             20.0 / 3.0 * fun[3] - 15.0 / 4.0 * fun[4] + 6.0 / 5.0 * fun[5] -
             1.0 / 6.0 * fun[6]);

  dfun[1] = (-49.0 / 20.0 * fun[1] + 6.0 * fun[2] - 15.0 / 2.0 * fun[3] +
             20.0 / 3.0 * fun[4] - 15.0 / 4.0 * fun[5] + 6.0 / 5.0 * fun[6] -
             1.0 / 6.0 * fun[7]);

  // 4-th order central difference
  for (int ir = 2; ir < endpoint - 2; ir++) {
    dfun[ir] =
        (fun[ir - 2] - 8.0 * fun[ir - 1] + 8.0 * fun[ir + 1] - fun[ir + 2]) /
        12.0;
  }

  // 6-th order backward difference for the last two points
  dfun[endpoint - 2] =
      -(-49.0 / 20.0 * fun[endpoint - 2] + 6.0 * fun[endpoint - 3] -
        15.0 / 2.0 * fun[endpoint - 4] + 20.0 / 3.0 * fun[endpoint - 5] -
        15.0 / 4.0 * fun[endpoint - 6] + 6.0 / 5.0 * fun[endpoint - 7] -
        1.0 / 6.0 * fun[endpoint - 8]);

  dfun[endpoint - 1] =
      -(-49.0 / 20.0 * fun[endpoint - 1] + 6.0 * fun[endpoint - 2] -
        15.0 / 2.0 * fun[endpoint - 3] + 20.0 / 3.0 * fun[endpoint - 4] -
        15.0 / 4.0 * fun[endpoint - 5] + 6.0 / 5.0 * fun[endpoint - 6] -
        1.0 / 6.0 * fun[endpoint - 7]);

  return dfun;
}

}  // namespace lsms

#endif  // LSMS_DIFF_HPP
