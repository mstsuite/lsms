/**
 *
 * RadialDFTLib
 *
 * RadialSrelSolver.cpp
 *
 * Created by Franco Moitzi on 3/3/23.
 *
 * Copyright (c) 2023 University of Leoben. All rights reserved.
 *
 */

#include "RadialSrelSolver.hpp"

void lsms::RadialSrelOut(std::complex<double> E, int kappa,
                           std::complex<double> *P, std::complex<double> *Q,
                           const double *R, const double *Rp, const double *V,
                           int end, int sign_correction, double *beta) {

  std::complex<double> yp[2];
  std::complex<double> f0[2];
  std::complex<double> f1[2];
  std::complex<double> f2[2];
  std::complex<double> f3[2];
  std::complex<double> delta;
  double lam;
  std::complex<double> i0;
  std::complex<double> i1;
  std::complex<double> M11;
  std::complex<double> M21;
  std::complex<double> M12;
  std::complex<double> M22;

  std::complex<double> k = kappa * kappa + kappa;

  double Z = -V[0] * R[0];

  if (Z > 0.0001) {
    *beta = std::sqrt(1 - (Z * Z) / (c * c) + (kappa * kappa + kappa));
  } else {
    *beta = abs(kappa);
  }

  double start_value = 1e-5;

  if (Z > 0.0001) {
    P[0] = start_value;
    Q[0] = start_value * (k)*c / (Z * (*beta + 1));
  } else {
    int l = -kappa - 1;
    P[0] = pow(start_value, (l + 1.0)) * pow(-1, sign_correction);
    Q[0] = (l + 1.0) * pow(start_value, l) * pow(-1, sign_correction);
  }

  // Predictor Corrector Euler
  cdirac_radial_func_nosoc(R[0], Rp[0], E, V[0], kappa, P[0], Q[0], f0);
  yp[0] = P[0] + f0[0];
  yp[1] = Q[0] + f0[1];
  cdirac_radial_func_nosoc(R[1], Rp[1], E, V[1], kappa, yp[0], yp[1], f1);
  P[1] = P[0] + (f0[0] + f1[0]) / 2.0;
  Q[1] = Q[0] + (f0[1] + f1[1]) / 2.0;

  // Adams Moulton 3
  cdirac_radial_func_nosoc(R[1], Rp[1], E, V[1], kappa, P[1], Q[1], f1);
  i0 = P[1] + (-f0[0] + 8.0 * f1[0]) / 12.0;
  i1 = Q[1] + (-f0[1] + 8.0 * f1[1]) / 12.0;

  lam = 5.0 / 12.0;
  delta = 1.0 + lam * lam * Rp[2] * Rp[2] *
                    ((E - V[2]) * (E - V[2]) / (c_2) + 2.0 * E - 2 * V[2] -
                     (k + 1.0) / (R[2] * R[2]));

  M11 = 1 + lam * Rp[2] / R[2];
  M12 = lam * Rp[2] * (E / c - V[2] / c + 2 * c);
  M21 = lam * Rp[2] *
        (-E / c + V[2] / c + (k) / (R[2] * R[2] * (E / c - V[2] / c + 2 * c)));
  M22 = 1 - lam * Rp[2] / R[2];
  P[2] = (M11 * i0 + M12 * i1) / delta;
  Q[2] = (M21 * i0 + M22 * i1) / delta;

  // Adams Moulton 4
  cdirac_radial_func_nosoc(R[2], Rp[2], E, V[2], kappa, P[2], Q[2], f2);
  i0 = P[2] + (f0[0] - 5.0 * f1[0] + 19.0 * f2[0]) / 24.0;
  i1 = Q[2] + (f0[1] - 5.0 * f1[1] + 19.0 * f2[1]) / 24.0;

  lam = 9.0 / 24.0;

  delta = 1.0 + lam * lam * Rp[3] * Rp[3] *
                    ((E - V[3]) * (E - V[3]) / (c_2) + 2.0 * E - 2 * V[3] -
                     (k + 1.0) / (R[3] * R[3]));

  M11 = 1 + lam * Rp[3] / R[3];
  M12 = lam * Rp[3] * (E / c - V[3] / c + 2 * c);
  M21 = lam * Rp[3] *
        (-E / c + V[3] / c + (k) / (R[3] * R[3] * (E / c - V[3] / c + 2 * c)));

  M22 = 1 - lam * Rp[3] / R[3];
  P[3] = (M11 * i0 + M12 * i1) / delta;
  Q[3] = (M21 * i0 + M22 * i1) / delta;

  // Adams Moulton 5
  for (int idx = 3; idx < end - 1; idx++) {
    cdirac_radial_func_nosoc(R[idx - 3], Rp[idx - 3], E, V[idx - 3], kappa,
                             P[idx - 3], Q[idx - 3], f0);

    cdirac_radial_func_nosoc(R[idx - 2], Rp[idx - 2], E, V[idx - 2], kappa,
                             P[idx - 2], Q[idx - 2], f1);

    cdirac_radial_func_nosoc(R[idx - 1], Rp[idx - 1], E, V[idx - 1], kappa,
                             P[idx - 1], Q[idx - 1], f2);

    cdirac_radial_func_nosoc(R[idx], Rp[idx], E, V[idx], kappa, P[idx], Q[idx],
                             f3);

    i0 = P[idx] +
         1.0 / 720 *
             (-19.0 * f0[0] + 106.0 * f1[0] - 264.0 * f2[0] + 646.0 * f3[0]);

    i1 = Q[idx] +
         1.0 / 720 *
             (-19.0 * f0[1] + 106.0 * f1[1] - 264.0 * f2[1] + 646.0 * f3[1]);

    lam = 251.0 / 720.0;

    delta = 1.0 + lam * lam * Rp[idx + 1] * Rp[idx + 1] *
                    ((E - V[idx + 1]) * (E - V[idx + 1]) / (c_2) + 2.0 * E -
                     2 * V[idx + 1] - (k + 1.0) / (R[idx + 1] * R[idx + 1]));

    M11 = 1 + lam * Rp[idx + 1] / R[idx + 1];

    M12 = lam * Rp[idx + 1] * (E / c - V[idx + 1] / c + 2 * c);

    M21 = lam * Rp[idx + 1] *
          (-E / c + V[idx + 1] / c +
           (k) / (R[idx + 1] * R[idx + 1] *
                  (E / c - V[idx + 1] / c +

                   2 * c)));

    M22 = 1 - lam * Rp[idx + 1] / R[idx + 1];

    P[idx + 1] = (M11 * i0 + M12 * i1) / delta;

    Q[idx + 1] = (M21 * i0 + M22 * i1) / delta;
  }
}

void lsms::RadialSrelIn(std::complex<double> E, int kappa,
                        std::complex<double> *P,
                        std::complex<double> *Q,
                        const double *R,
                        const double *Rp,
                        const double *V,
                        int end,
                        std::complex<double> *P_last,
                        std::complex<double> *Q_last,
                        int stop) {

  int idx = end - 1;
  int start_imax = idx;

  std::complex<double>  f0[2], f1[2], f2[2];

  std::complex<double>  delta, lam, i0, i1;

  std::complex<double>  M11;
  std::complex<double>  M21;
  std::complex<double>  M12;
  std::complex<double>  M22;

  const double min_tol = 1.0e-10;

  std::complex<double>  lambda = std::sqrt(-2.0 * E - E * E / (c_2));

  // Use asymptoics for starting
  for (idx = end - 1; idx >= stop; idx--) {
    P[idx] = exp(-lambda * R[idx]) / sqrt(-E / (E + 2 * c_2));
    Q[idx] = -exp(-lambda * R[idx]);

    if (std::abs(P[idx]) > min_tol) {
      start_imax = idx;
      break;
    }
  }

  /*
   * Calculate 3 more points and then use predictor corrector to improve the
   * last point
   */

  P[start_imax - 1] =
      exp(-lambda * R[start_imax - 1]) / sqrt(-E / (E + 2 * c_2));

  Q[start_imax - 1] = -exp(-lambda * R[start_imax - 1]);

  P[start_imax - 2] =
      exp(-lambda * R[start_imax - 2]) / sqrt(-E / (E + 2 * c_2));

  Q[start_imax - 2] = -exp(-lambda * R[start_imax - 2]);

  for (idx = start_imax; idx > stop + 3; --idx) {

    cdirac_radial_func_nosoc(R[idx], Rp[idx], E, V[idx], kappa,
                             P[idx], Q[idx], f0);

    cdirac_radial_func_nosoc(R[idx - 1], Rp[idx - 1], E, V[idx - 1], kappa,
                             P[idx - 1], Q[idx - 1], f1);

    cdirac_radial_func_nosoc(R[idx - 2], Rp[idx - 2], E, V[idx - 2], kappa,
                             P[idx - 2], Q[idx - 2], f2);

    i0 = P[idx - 2] - 1.0 / 24 * (f0[0] - 5.0 * f1[0] + 19.0 * f2[0]);
    i1 = Q[idx - 2] - 1.0 / 24 * (f0[1] - 5.0 * f1[1] + 19.0 * f2[1]);
    lam = -9.0 / 24.0;

    delta = 1.0 +
        lam * lam * Rp[idx - 3] * Rp[idx - 3] *
            (-kappa * kappa / (R[idx - 3] * R[idx - 3]) -
                (-E + V[idx - 3]) * (2 * c + (E - V[idx - 3]) / c) / c);

    //M11 = -kappa * lam * Rp[idx - 3] / R[idx - 3] + 1;
    //M22 = kappa * lam * Rp[idx - 3] / R[idx - 3] + 1;
    //M21 = lam * Rp[idx - 3] * (-E + V[idx - 3]) / c;
    //M12 = lam * Rp[idx - 3] * (2 * c + (E - V[idx - 3]) / c);

    P[idx - 3] = (M11 * i0 + M12 * i1) / delta;
    Q[idx - 3] = (M21 * i0 + M22 * i1) / delta;
  }

  idx = stop + 3;

  cdirac_radial_func_nosoc(R[idx], Rp[idx], E, V[idx], kappa,
                           P[idx], Q[idx], f0);

  cdirac_radial_func_nosoc(R[idx - 1], Rp[idx - 1], E, V[idx - 1], kappa,
                           P[idx - 1], Q[idx - 1], f1);

  cdirac_radial_func_nosoc(R[idx - 2], Rp[idx - 2], E, V[idx - 2], kappa,
                           P[idx - 2], Q[idx - 2], f2);

  i0 = P[idx - 2] - 1.0 / 24 * (f0[0] - 5.0 * f1[0] + 19.0 * f2[0]);

  i1 = Q[idx - 2] - 1.0 / 24 * (f0[1] - 5.0 * f1[1] + 19.0 * f2[1]);

  lam = -9.0 / 24.0;

//  delta = 1 +
//      lam * lam * Rp[idx - 3] * Rp[idx - 3] *
//          (-kappa * kappa / (R[idx - 3] * R[idx - 3]) -
//              (-E + V[idx - 3]) * (2 * c + (E - V[idx - 3]) / c) / c);

//  M11 = -kappa * lam * Rp[idx - 3] / R[idx - 3] + 1;
//
//  M22 = kappa * lam * Rp[idx - 3] / R[idx - 3] + 1;
//
//  M21 = lam * Rp[idx - 3] * (-E + V[idx - 3]) / c;
//
//  M12 = lam * Rp[idx - 3] * (2 * c + (E - V[idx - 3]) / c);

  *P_last = (M11 * i0 + M12 * i1) / delta;

  *Q_last = (M21 * i0 + M22 * i1) / delta;
}
