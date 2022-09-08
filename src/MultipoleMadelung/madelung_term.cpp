//
// Created by F.Moitzi on 07.09.2022.
//

#include "madelung_term.hpp"

#include "Indices.hpp"
#include "Coeficients.hpp"

#include "spherical_harmonics.hpp"

//static void spherical_harmonics(std::vector<Real> &vec, int lmax,
//                                std::vector<Complex> &ylm,
//                                std::vector<Real> &plm) {
//  using namespace std::complex_literals;
//
//  auto q2 = vec[0] * vec[0] + vec[1] * vec[1];
//  auto r = sqrt(q2 + vec[2] * vec[2]);
//  auto q = sqrt(q2);
//
//  Complex iphi = 1i * atan2(vec[1] / q, vec[0] / q);
//
//  associatedLegendreFunctionNormalized(vec[2] / r, lmax, plm.data());
//
//  int kl = 0;
//
//  for (int l = 0; l <= lmax; l++) {
//    ylm[kl] = plm[plmIdx(l, 0)];
//
//    for (int m = 1; m <= l; m++) {
//      ylm[kl + 1] = plm[plmIdx(l, m)] * std::exp(iphi * (Complex)m);
//      ylm[kl - 1] = std::conj(ylm[kl + 1]) * std::pow(-1, m);
//    }
//
//    kl += (l + 1) * 2;
//  }
//}

void lsms::dlsum(std::vector<Real> &aij, matrix<Real> &rslat, int nrslat,
                 int ibegin, matrix<Real> &knlat, int nknlat, double omega,
                 int lmax_mad, double eta, std::vector<Complex> &dlm) {
  auto kmax_mad = (lmax_mad + 1) * (lmax_mad + 1);
  std::vector<Complex> Ylm((lmax_mad + 1) * (lmax_mad + 1), Complex(0.0, 0.0));
  std::vector<Real> Plm((lmax_mad + 1) * (lmax_mad + 2) / 2, 0.0);

  std::vector<Real> vec(3);

  std::fill(dlm.begin(), dlm.end(), std::complex<double>(0,0));

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

    for (auto kl = kmax_mad - 1; kl > 0; kl--) {
      auto l = AngularMomentumIndices::lofk[kl];
      dlm[kl] = dlm[kl] + gamma_l[l] * Ylm[kl] / std::pow(vhalf, l + 1);
    }
  }

  auto rfac = 4.0 * std::sqrt(M_PI);
  for (int i = 1; i < kmax_mad; ++i) {
    dlm[i] *= rfac;
  }

  rfac = -eta * eta / 4.0;

  std::vector<Complex> ctmp(kmax_mad, 0.0);

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

    Real tfac, sintfac, costfac;

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

}
