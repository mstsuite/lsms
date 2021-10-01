//
// Created by F.Moitzi on 28.09.2021.
//

#ifndef LSMS_POISSON_HPP
#define LSMS_POISSON_HPP


#define _USE_MATH_DEFINES

#include <limits>
#include <cmath>
#include <vector>

#include "integrator.hpp"

namespace lsms {


  template<typename T>
  T static simps_r(const T *val,
                   const T *R,
                   const T *Rp,
                   const int end) {

    T intval = 0.0;
    int idx;

    for (idx = 0; idx < end - 2; idx += 2) {
      intval += 1.0 / 3.0 * (val[idx] * R[idx] * Rp[idx] +
                             4 * val[idx + 1] * R[idx + 1] * Rp[idx + 1] +
                             val[idx + 2] * R[idx + 2] * Rp[idx + 2]);
    };

    if ((end % 2) == 0) {
      intval += 0.5 * (val[idx] * R[idx] * Rp[idx] + val[idx + 1] * R[idx + 1] * Rp[idx + 1]);
    };
    return intval;
  }

  template<typename T>
  static inline void rad_func(
      const double &r,
      const double &rp,
      const T &n,
      const T &q,
      T(&f)[2]) {

    f[0] = rp * q;
    f[1] = rp * (-n / (r * r) - 2.0 / r * q);
  }

  /**
   * Radial poisson solver
   *
   * Density will be rho(r) = 4 * PI * r * r * n(r)
   *
   */
  template<typename T>
  void radial_poisson(std::vector<T> &VHartree,
                      std::vector<T> &VHartreeDeriv,
                      const std::vector<T> &R,
                      const std::vector<T> &Rp,
                      const std::vector<T> &density,
                      int end) {

    int idx;
    T yp[2];
    T f0[2];
    T f1[2];
    T f2[2];
    T f3[2];
    T f4[2];
    T delta;
    T lam;
    T i0;
    T i1;


    std::vector<T> density_r(end, 0.0);
    for (int i = 0; i < end; i++) {
      density_r[i] = density[i] / R[i];
    }

    // Starter condition
    VHartree[0] = radialIntegral(density_r, R, end);
    VHartreeDeriv[0] = 0;


    // PC Euler
    rad_func(R[0], Rp[0], density[0], VHartreeDeriv[0], f0);
    yp[0] = VHartree[0] + f0[0];
    yp[1] = VHartreeDeriv[0] + f0[1];
    rad_func(R[1], Rp[1], density[1], yp[1], f1);
    VHartree[1] = VHartree[0] + (f0[0] + f1[0]) / 2;
    VHartreeDeriv[1] = VHartreeDeriv[0] + (f0[1] + f1[1]) / 2;

  /*
    for (int idx = 1; idx < end; ++idx) {
      rad_func(R[idx - 1], Rp[idx - 1], density[idx - 1], VHartreeDeriv[idx - 1], f0);
      yp[0] = VHartree[idx - 1] + f0[0];
      yp[1] = VHartreeDeriv[idx - 1] + f0[1];
      rad_func(R[idx], Rp[idx], density[idx], yp[1], f1);
      VHartree[idx] = VHartree[idx - 1] + (f0[0] + f1[0]) / 2;
      VHartreeDeriv[idx] = VHartreeDeriv[idx - 1] + (f0[1] + f1[1]) / 2;

    }
    return;
*/

    // Adams Moulton 3
    idx = 1;

    rad_func(R[idx], Rp[idx],
             density[idx], VHartreeDeriv[idx], f1);
    lam = 5.0 / 12.0;

    i0 = VHartree[idx] + 1.0 / 12.0 * (-f0[0] + 8.0 * f1[0]);
    i1 = VHartreeDeriv[idx] + 1.0 / 12.0 * (-f0[1] + 8.0 * f1[1])
         - Rp[idx + 1] * lam * density[idx + 1] / (R[idx + 1] * R[idx + 1]);

    delta = 1 + lam * 2.0 * Rp[idx + 1] / R[idx + 1];

    VHartree[idx + 1] = i0 + lam * Rp[idx + 1] * i1 / delta;
    VHartreeDeriv[idx + 1] = i1 / delta;

    // Adams Moulton 4
    idx = 2;

    rad_func(R[idx], Rp[idx],
             density[idx], VHartreeDeriv[idx], f2);

    lam = 9.0 / 24.0;

    i0 = VHartree[idx] + 1.0 / 24.0 * (f0[0] - 5.0 * f1[0] + 19.0 * f2[0]);
    i1 = VHartreeDeriv[idx] + 1.0 / 24.0 * (f0[1] - 5.0 * f1[1] + 19.0 * f2[1])
         - Rp[idx + 1] * lam * density[idx + 1] / (R[idx + 1] * R[idx + 1]);

    delta = 1 + lam * 2.0 * Rp[idx + 1] / R[idx + 1];

    VHartree[idx + 1] = i0 + lam * Rp[idx + 1] * i1 / delta;
    VHartreeDeriv[idx + 1] = i1 / delta;

    // Adams Moulton 5
    idx = 3;
    rad_func(R[idx], Rp[idx],
             density[idx], VHartreeDeriv[idx], f3);

    lam = 251.0 / 720.0;

    i0 = VHartree[idx] + 1.0 / 720.0 * (-19.0 * f0[0] - 106.0 * f1[0] - 264 * f2[0] + 646 * f3[0]);
    i1 = VHartreeDeriv[idx] + 1.0 / 720.0 * (-19.0 * f0[1] - 106.0 * f1[1] - 264.0 * f2[0] + 646 * f3[1])
         - Rp[idx + 1] * lam * density[idx + 1] / (R[idx + 1] * R[idx + 1]);

    delta = 1 + lam * 2.0 * Rp[idx + 1] / R[idx + 1];

    VHartree[idx + 1] = i0 + lam * Rp[idx + 1] * i1 / delta;
    VHartreeDeriv[idx + 1] = i1 / delta;

    // Adams Moulton 6
    for (idx = 4; idx < end - 1; ++idx) {

      rad_func(R[idx - 4], Rp[idx - 4],
               density[idx - 4], VHartreeDeriv[idx - 4], f0);
      rad_func(R[idx - 3], Rp[idx - 3],
               density[idx - 3], VHartreeDeriv[idx - 3], f1);
      rad_func(R[idx - 2], Rp[idx - 2],
               density[idx - 2], VHartreeDeriv[idx - 2], f2);
      rad_func(R[idx - 1], Rp[idx - 1],
               density[idx - 1], VHartreeDeriv[idx - 1], f3);
      rad_func(R[idx], Rp[idx],
               density[idx], VHartreeDeriv[idx], f4);

      lam = 475.0 / 1440.0;

      i0 = VHartree[idx] + 1.0 / 1440.0 * (27.0 * f0[0] - 173.0 * f1[0]
                                           + 482.0 * f2[0] - 798.0 * f3[0] + 1427.0 * f4[0]);

      i1 = VHartreeDeriv[idx] + 1.0 / 1440.0 * (27.0 * f0[1] - 173.0 * f1[1]
                                                + 482.0 * f2[1] - 798.0 * f3[1] + 1427.0 * f4[1])
           - Rp[idx + 1] * lam * density[idx + 1] / (R[idx + 1] * R[idx + 1]);

      delta = 1 + lam * 2.0 * Rp[idx + 1] / R[idx + 1];

      VHartree[idx + 1] = i0 + lam * Rp[idx + 1] * i1 / delta;
      VHartreeDeriv[idx + 1] = i1 / delta;

    }

  }
}


#endif //LSMS_POISSON_HPP
