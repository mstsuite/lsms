
#ifndef LSMS_RADIAL_SOLVER_H
#define LSMS_RADIAL_SOLVER_H

#include <cfloat>
#include <cstdlib>
#include <cmath>

namespace lsms {

constexpr int ERROR_WRONG_NUMBER_NODES = -3;
constexpr int ERROR_STILL_BISECTION = -2;
constexpr int ERROR_NOT_REACHED_ACCURACY = -1;

constexpr double U_LIM = 1e5;
constexpr double L_LIM = 1e-5;
constexpr double C = 137.0359991; // Speed of light in atomic units (inverse \f$\alpha\f$)
constexpr double C_SQ = 18778.86504933520081;

template<typename T, typename D>
inline void rel_rad_func(T r, T rp,
                         T E, T vpot,
                         D k_qn, T p,
                         T q, T val[2]) {
  val[0] = (-k_qn / r * p + ((E - vpot) / C + 2 * C) * q) * rp;
  val[1] = (-((E - vpot) / C) * p + k_qn / r * q) * rp;
}

template<typename T, typename D>
inline void nrel_rad_func(T r, T rp,
                          T E, T vpot,
                          D l, T p, T q,
                          T val[2]) {
  val[0] = q * rp;
  val[1] = 2.0 * ((vpot - E + (l * l + l) / (2.0 * r * r)) * p) * rp;
}

template<typename T, typename D, typename V>
inline V calculate_ctp(T E, D l, const T *R,
                       const T *pot, V end) {
  auto idx = end;

  for (idx = end - 1; idx >= 0; idx--) {
    if ((E - pot[idx] - (l * l + l) / (2.0 * R[idx] * R[idx])) > 0) {
      return idx;
    }
  }

  return end;
}

template<typename T, typename D, typename V>
inline int count_sign_sign_switch(const T *P, D start,
                                  V end) {
  auto count = 0;

  for (auto idx = start; idx < end - 1; idx++) {
    if (P[idx] * P[idx + 1] < 0) {
      count += 1;
    }
  }

  return count;
}

template<typename T, typename D>
inline D last_valid_point(const T *P, D end) {
  auto i = end - 1;

  for (i = 0; i < end; i++) {
    if (std::fabs(P[i]) >= lsms::U_LIM) {
      break;
    }
  }

  return i;

}

void nrel_out(double Z, double energy, int l_qn, double *p_array, double *q_array,
              const double *__restrict__ r_mesh, double h,
              const double *__restrict__ pot, std::size_t end, int sign_switch);

void nrel_in(double energy, int l_qn, double *p_array, double *q_array,
             const double *__restrict__ r_mesh, double h,
             const double *__restrict__ pot, std::size_t end, double &p_last,
             double &q_last, std::size_t stop, std::size_t &imax);

void rel_in(double energy, int k_qn, double *p_array, double *q_array,
            const double *__restrict__ r_mesh, double h,
            const double *__restrict__ pot, std::size_t end, double &p_last,
            double &q_last, std::size_t stop, std::size_t &imax);

void rel_out(double energy, int k_qn, double *p_array, double *q_array,
             const double *__restrict__ r_mesh, double h,
             const double *__restrict__ pot, std::size_t end, int sign_switch);

inline double rel_integration(const double *__restrict__ p_array,
                              const double *__restrict__ q_array,
                              const double *__restrict__ r_mesh,
                              double h, std::size_t end);

inline double nrel_integration(const double *__restrict__ P,
                               const double *__restrict__ r_mesh,
                               double h,
                               std::size_t end);

double nonrel_eigenenergies(double Z, int n_qn, int l_qn, int k_qn, double *dens,
                            double *p_array, double *q_array, const double *__restrict__ r_mesh,
                            double h,
                            const double *__restrict__ pot, std::size_t end,
                            double energy_rel_tol, int max_iter, double e_init,
                            int &converged, double &delta_energy);

double rel_eigenenergies(
    double Z, int n_qn, int l_qn, int k_qn, double *dens, double *p_array, double *q_array,
    const double *__restrict__ r_mesh, double h,
    const double *__restrict__ pot, std::size_t end, double energy_rel_tol, int max_iter,
    double e_init, int &converged, double &delta_energy);

double rel_energy_start(int n, int k_qn, double Z);

double nrel_energy_start(int n, double Z);

}  // namespace lsms

#endif  // LSMS_RADIAL_SOLVER_H
