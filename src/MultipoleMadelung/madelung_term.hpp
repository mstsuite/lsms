//
// Created by F.Moitzi on 07.01.2022.
//

#ifndef SRC_MADELUNG_MADELUNG_TERM_HPP
#define SRC_MADELUNG_MADELUNG_TERM_HPP

#include <algorithm>
#include <complex>
#include <numeric>
#include <vector>

#include "common.hpp"
#include "utils.hpp"

namespace lsms {

/**
 * Reciprocal space term of Madelung sum:
 *
 *                                       2   2       -> ->
 *                4*pi          1    -eta *Kq /4 - i*Kq*aij
 *       term1 =  ---- * sum  ----- e
 *                tau    q<>0    2
 *                             Kq
 *
 *       note sum starts at 2 since kn=0.0 of 1/(rs-aij) is
 *       canceled by kn=0.0 of the 1/rs sum.
 *
 */
  template<class T>
  inline T reciprocal_space_term(matrix<T> &knlat, std::vector<T> &knlatsq,
                                 std::vector<T> &aij, int nknlat, double eta,
                                 double omega) {
    auto term12 = 0.0;
    auto fac = eta * eta / 4.0;

    for (auto i = nknlat - 1; i > 0; i--) {
      term12 += std::exp(-fac * knlatsq[i]) / knlatsq[i] *
                std::cos(knlat(0, i) * aij[0] + knlat(1, i) * aij[1] +
                         knlat(2, i) * aij[2]);
    }
    term12 = 4.0 * M_PI / omega * term12;

    return term12;
  }

/**
 * Real space term of Madelung sum:
 *
 *                              ->   ->
 *                     1 - erf(|Rn + aij|/eta)
 *        term2 = sum -------------------------
 *                 n         ->   ->
 *                         |Rn + aij|
 *
 *        note for calculation of aij=0.0 term ibegin=2.
 *
 */
  template<class T>
  inline T real_space_term(matrix<T> &rslat, std::vector<T> &aij, int nrslat,
                           int ibegin, T eta) {
    /*
     *  subtract aij from rslat and calculate rslatmd which is used in
     *  calculating the real-space integral
     *  rslatmd, and aij are in the units of a0 = 1
     */

    auto rterm = 0.0;

    double rslatmd;

    for (auto i = ibegin; i < nrslat; i++) {
      rslatmd = std::sqrt((rslat(0, i) - aij[0]) * (rslat(0, i) - aij[0]) +
                          (rslat(1, i) - aij[1]) * (rslat(1, i) - aij[1]) +
                          (rslat(2, i) - aij[2]) * (rslat(2, i) - aij[2]));

      rterm += std::erfc(rslatmd / eta) / rslatmd;
    }

    return rterm;
  }


/**
 *
 * Dl sum
 *
 */
  template<class T>
  std::vector<std::complex<T>> dlsum(std::vector<T> &aij,
                                     matrix<T> &rslat,
                                     int nrslat,
                                     int ibegin,
                                     matrix<T> &knlat,
                                     int nknlat,
                                     double omega,
                                     int lmax_mad,
                                     int kmax_mad, double eta) {

    std::vector<std::complex<T>> Ylm(kmax_mad, std::complex<T>(0.0, 0.0));
    std::vector<T> vec(3);
    std::vector<std::complex<T>> dlm(kmax_mad, 0.0);

    auto aij2 = std::inner_product(aij.begin(), aij.end(), aij.begin(), 0.0);

    for (int i = nrslat - 1; i >= ibegin; i--) {
      for (int j = 0; j < 3; j++) {
        vec[j] = aij[j] + rslat(j, i);
      }

      auto vlen = norm(vec.begin(), vec.end());

      // Ylm
      sph_harm_1(vec.data(), &lmax_mad, Ylm.data());

      // Gamma
      auto gamma_l = gamma_func(vlen / eta, lmax_mad);

      auto vhalf = 0.5 * vlen;


      auto kl = 0;

      for (auto l = 0; l <= lmax_mad; l++) {
        for (auto m = -l; m <= l; m++) {
          dlm[kl] = dlm[kl] + gamma_l[l] * Ylm[kl] / std::pow(vhalf, l + 1);
          kl++;
        }
      }

    }

    auto rfac = 4.0 * std::sqrt(M_PI);
    for (int i = 1; i < kmax_mad; ++i) {
      dlm[i] *= rfac;
    }

    rfac = -eta * eta / 4.0;

    std::vector<std::complex<T>> ctmp(kmax_mad, 0.0);

    for (int i = nknlat - 1; i > 0; i--) {
      for (int j = 0; j < 3; j++) {
        vec[j] = knlat(j, i);
      }

      auto vlen = norm(vec.begin(), vec.end());
      auto knlatsq = norm_sq(vec.begin(), vec.end());

      // std::printf("%16.12f %16.12f %16.12f %16.12f %16.12f\n", vec[0], vec[1],
      // vec[2], vlen, knlatsq);

      // Ylm
      sph_harm_1(vec.data(), &lmax_mad, Ylm.data());

      auto expfac = std::exp(rfac * knlatsq) / knlatsq;

      T tfac;
      T sintfac, costfac;

      if (aij2 > 1e-8) {
        tfac =
            -(knlat(0, i) * aij[0] + knlat(1, i) * aij[1] + knlat(2, i) * aij[2]);
        sintfac = sin(tfac);
        costfac = cos(tfac);
      } else {
        tfac = 0.0;
        sintfac = 0.0;
        costfac = 1.0;
      }

      auto cfac = -vlen;
      for (auto l = 1; l <= lmax_mad; l += 2) {
        for (auto m = -l; m <= l; m++) {
          auto kl = (l + 1) * (l + 1) - l + m - 1;
          ctmp[kl] = ctmp[kl] + expfac * cfac * sintfac * Ylm[kl];
        }
        cfac = -cfac * knlatsq;
      }

      cfac = -knlatsq;
      for (auto l = 2; l <= lmax_mad; l += 2) {
        for (auto m = -l; m <= l; m++) {
          auto kl = (l + 1) * (l + 1) - l + m - 1;
          ctmp[kl] = ctmp[kl] + expfac * cfac * costfac * Ylm[kl];
        }
        cfac = -cfac * knlatsq;
      }
    }

    rfac = 16.0 * M_PI * M_PI / omega;

    for (auto i = 1; i < kmax_mad; i++) {
      dlm[i] = dlm[i] + rfac * ctmp[i];
    }

    return dlm;
  }


}  // namespace lsms

#endif  // SRC_MADELUNG_MADELUNG_TERM_HPP
