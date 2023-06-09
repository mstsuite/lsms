/**
 *
 *
 * RadialSrelSolver.hpp
 *
 * Created by Franco Moitzi on 3/3/23.
 *
 */
#ifndef LSMS_RADIALSRELSOLVER_HPP
#define LSMS_RADIALSRELSOLVER_HPP

#include <complex>

namespace lsms {

constexpr double c = 137.0359991;
constexpr double c_2 = 18778.86504933520081;

inline void cdirac_radial_func_nosoc(const double r, const double rp,
                                     const std::complex<double> E,
                                     const double vpot, const int kappa,
                                     const std::complex<double> p,
                                     const std::complex<double> q,
                                     std::complex<double> f[2]) {

  std::complex<double> k = kappa * kappa + kappa;

  f[0] = rp * (p / r + ((E - vpot) / c + 2 * c) * q);
  f[1] =
      rp * (-q / r +
          (-E / c + vpot / c + k / (r * r * ((E - vpot) / c + 2 * c))) * p);
}

void RadialSrelOut(std::complex<double> E, int kappa, std::complex<double> *P,
                   std::complex<double> *Q, const double *R, const double *Rp,
                   const double *V, int end, int sign_correction,
                   double *beta);

void RadialSrelIn(std::complex<double> E, int kappa,
                  std::complex<double> *P,
                  std::complex<double> *Q,
                  const double *R,
                  const double *Rp,
                  const double *V,
                  int end,
                  std::complex<double> *P_last,
                  std::complex<double> *Q_last,
                  int stop);

} // namespace lsms

#endif // LSMS_RADIALSRELSOLVER_HPP
