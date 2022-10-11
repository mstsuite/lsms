

#include "radialSolver.hpp"

#include <cfloat>
#include <cstdlib>
#include <cmath>

void lsms::nrel_in(double energy, int l_qn, double *p_array, double *q_array,
                   const double *__restrict__ r_mesh,
                   double h,
                   const double *__restrict__ pot, std::size_t end, double &p_last,
                   double &q_last, std::size_t stop, std::size_t &imax) {
  auto idx = end - 1;
  auto start_idx = idx;
  imax = start_idx;

  double val0[2], val1[2], val2[2], val3[2];

  double delta;
  double chi, inv0, inv1, step_val1, step_val3, step_val2;

  constexpr auto tolerance = DBL_MIN;
  double omega = std::sqrt(-2 * energy);

  for (idx = end - 1; idx >= stop; idx--) {
    p_array[idx] = std::exp(-omega * r_mesh[idx]);
    q_array[idx] = -std::exp(-omega * r_mesh[idx]) * omega;

    if (std::fabs(p_array[idx]) > tolerance && std::fabs(q_array[idx]) > tolerance) {
      start_idx = idx;
      break;
    } else {
      p_array[idx] = 0.0;
      q_array[idx] = 0.0;
    }
  }

  imax = start_idx;

  p_array[start_idx - 1] = std::exp(-omega * r_mesh[start_idx - 1]);
  q_array[start_idx - 1] = -std::exp(-omega * r_mesh[start_idx - 1]) * omega;
  p_array[start_idx - 2] = std::exp(-omega * r_mesh[start_idx - 2]);
  q_array[start_idx - 2] = -std::exp(-omega * r_mesh[start_idx - 2]) * omega;

  idx = start_idx;
  nrel_rad_func(r_mesh[idx], h * r_mesh[idx], energy, pot[idx], l_qn, p_array[idx], q_array[idx], val0);
  nrel_rad_func(r_mesh[idx - 1], h * r_mesh[idx - 1], energy, pot[idx - 1], l_qn, p_array[idx - 1],
                q_array[idx - 1], val1);
  nrel_rad_func(r_mesh[idx - 2], h * r_mesh[idx - 2], energy, pot[idx - 2], l_qn, p_array[idx - 2],
                q_array[idx - 2], val2);

  inv0 = p_array[idx - 2] - 1.0 / 24 * (val0[0] - 5.0 * val1[0] + 19.0 * val2[0]);
  inv1 = q_array[idx - 2] - 1.0 / 24 * (val0[1] - 5.0 * val1[1] + 19.0 * val2[1]);
  chi = -9.0 / 24.0;

  step_val1 = r_mesh[idx - 3];
  step_val2 = h * r_mesh[idx - 3];
  step_val3 = -2 * h * r_mesh[idx - 3] * (energy - pot[idx - 3] - l_qn * (l_qn + 1.0) / (2 * step_val1 * step_val1));
  delta = 1 - chi * chi * step_val2 * step_val3;

  p_array[idx - 3] = (inv0 + chi * step_val2 * inv1) / delta;
  q_array[idx - 3] = (chi * step_val3 * inv0 + inv1) / delta;

  
  for (idx = start_idx; idx > stop + 4; --idx) {
    nrel_rad_func(r_mesh[idx - 3], h * r_mesh[idx - 3], energy, pot[idx - 3], l_qn, p_array[idx - 3],
                  q_array[idx - 3], val3);

    inv0 = p_array[idx - 3] -
        1.0 / 720 *
            (-19.0 * val0[0] + 106.0 * val1[0] - 264.0 * val2[0] + 646.0 * val3[0]);

    inv1 = q_array[idx - 3] -
        1.0 / 720 *
            (-19.0 * val0[1] + 106.0 * val1[1] - 264.0 * val2[1] + 646.0 * val3[1]);

    chi = -251.0 / 720.0;

    step_val1 = r_mesh[idx - 4];
    step_val2 = h * r_mesh[idx - 4];

    step_val3 = -2 * h * r_mesh[idx - 4] * (energy - pot[idx - 4] - l_qn * (l_qn + 1.0) / (2 * step_val1 * step_val1));
    delta = 1 - chi * chi * step_val2 * step_val3;

    p_array[idx - 4] = (inv0 + chi * step_val2 * inv1) / delta;
    q_array[idx - 4] = (chi * step_val3 * inv0 + inv1) / delta;

    val0[0] = val1[0];
    val0[1] = val1[1];

    val1[0] = val2[0];
    val1[1] = val2[1];

    val2[0] = val3[0];
    val2[1] = val3[1];
  }

  idx = stop + 4;
  nrel_rad_func(r_mesh[idx], h * r_mesh[idx], energy, pot[idx], l_qn, p_array[idx], q_array[idx], val0);
  nrel_rad_func(r_mesh[idx - 1], h * r_mesh[idx - 1], energy, pot[idx - 1], l_qn, p_array[idx - 1],
                q_array[idx - 1], val1);
  nrel_rad_func(r_mesh[idx - 2], h * r_mesh[idx - 2], energy, pot[idx - 2], l_qn, p_array[idx - 2],
                q_array[idx - 2], val2);
  nrel_rad_func(r_mesh[idx - 3], h * r_mesh[idx - 3], energy, pot[idx - 3], l_qn, p_array[idx - 3],
                q_array[idx - 3], val3);

  inv0 = p_array[idx - 3] -
      1.0 / 720 *
          (-19.0 * val0[0] + 106.0 * val1[0] - 264.0 * val2[0] + 646.0 * val3[0]);

  inv1 = q_array[idx - 3] -
      1.0 / 720 *
          (-19.0 * val0[1] + 106.0 * val1[1] - 264.0 * val2[1] + 646.0 * val3[1]);

  chi = -251.0 / 720.0;

  step_val1 = r_mesh[idx - 4];
  step_val2 = h * r_mesh[idx - 4];

  step_val3 = -2 * h * r_mesh[idx - 4] * (energy - pot[idx - 4] - l_qn * (l_qn + 1.0) / (2 * step_val1 * step_val1));

  delta = 1 - chi * chi * step_val2 * step_val3;

  p_last = (inv0 + chi * step_val2 * inv1) / delta;
  q_last = (chi * step_val3 * inv0 + inv1) / delta;
}

void lsms::nrel_out(double Z, double energy, int l_qn, double *p_array, double *q_array,
                    const double *__restrict__ r_mesh,
                    double h,
                    const double *__restrict__ pot, std::size_t end,
                    int sign_switch) {
  auto idx = end - end;
  double yp[2], val0[2], val1[2], val2[2], val3[2];
  double delta, step_val3, chi, step_val1, step_val2, inv0, inv1;

  if (l_qn == 0) {
    if (Z <= L_LIM) {
      p_array[0] = (2 + 2 * energy * r_mesh[0] * r_mesh[0]) * std::pow(-1, sign_switch);
      q_array[0] = energy * r_mesh[0] * std::pow(-1, sign_switch);
    } else {
      p_array[0] = (r_mesh[0] - Z * r_mesh[0] * r_mesh[0]) * std::pow(-1, sign_switch);
      q_array[0] = (1 - 2 * r_mesh[0] * Z) * std::pow(-1, sign_switch);
    }
  } else {
    p_array[0] = std::pow(r_mesh[0], (l_qn + 1.0)) * std::pow(-1, sign_switch);
    q_array[0] = (l_qn + 1.0) * std::pow(r_mesh[0], l_qn) * std::pow(-1, sign_switch);
  }

  nrel_rad_func(r_mesh[0], h * r_mesh[0], energy, pot[0], l_qn, p_array[0], q_array[0], val0);
  yp[0] = p_array[0] + val0[0];
  yp[1] = q_array[0] + val0[1];

  nrel_rad_func(r_mesh[1], h * r_mesh[1], energy, pot[1], l_qn, yp[0], yp[1], val1);
  p_array[1] = p_array[0] + (val0[0] + val1[0]) / 2;
  q_array[1] = q_array[0] + (val0[1] + val1[1]) / 2;

  nrel_rad_func(r_mesh[1], h * r_mesh[1], energy, pot[1], l_qn, p_array[1], q_array[1], val1);

  inv0 = p_array[1] + 1.0 / 12.0 * (-val0[0] + 8.0 * val1[0]);
  inv1 = q_array[1] + 1.0 / 12.0 * (-val0[1] + 8.0 * val1[1]);

  chi = 5.0 / 12.0;
  step_val3 = -2 * h * r_mesh[2] * (energy - pot[2] - l_qn * (l_qn + 1.0) / (2 * r_mesh[2] * r_mesh[2]));

  delta = 1 - chi * chi * h * r_mesh[2] * step_val3;
  p_array[2] = (inv0 + chi * h * r_mesh[2] * inv1) / delta;
  q_array[2] = (chi * step_val3 * inv0 + inv1) / delta;

  nrel_rad_func(r_mesh[2], h * r_mesh[2], energy, pot[2], l_qn, p_array[2], q_array[2], val2);
  inv0 = p_array[2] + 1.0 / 24 * (val0[0] - 5.0 * val1[0] + 19.0 * val2[0]);
  inv1 = q_array[2] + 1.0 / 24 * (val0[1] - 5.0 * val1[1] + 19.0 * val2[1]);

  chi = 9.0 / 24.0;

  step_val1 = r_mesh[3];
  step_val2 = h * r_mesh[3];

  step_val3 = -2 * h * r_mesh[3] * (energy - pot[3] - l_qn * (l_qn + 1.0) / (2 * step_val1 * step_val1));
  delta = 1 - chi * chi * step_val2 * step_val3;

  p_array[3] = (inv0 + chi * step_val2 * inv1) / delta;
  q_array[3] = (chi * step_val3 * inv0 + inv1) / delta;

  
  for (idx = 3; idx < end - 1; ++idx) {
    nrel_rad_func(r_mesh[idx], h * r_mesh[idx], energy, pot[idx], l_qn, p_array[idx], q_array[idx], val3);

    inv0 = p_array[idx] +
        1.0 / 720 *
            (-19.0 * val0[0] + 106.0 * val1[0] - 264.0 * val2[0] + 646.0 * val3[0]);

    inv1 = q_array[idx] +
        1.0 / 720 *
            (-19.0 * val0[1] + 106.0 * val1[1] - 264.0 * val2[1] + 646.0 * val3[1]);

    chi = 251.0 / 720.0;

    step_val1 = r_mesh[idx + 1];
    step_val2 = h * r_mesh[idx + 1];

    step_val3 = -2 * h * r_mesh[idx + 1] * (energy - pot[idx + 1] - l_qn * (l_qn + 1.0) / (2 * step_val1 * step_val1));

    delta = 1 - chi * chi * step_val2 * step_val3;

    p_array[idx + 1] = (inv0 + chi * step_val2 * inv1) / delta;
    q_array[idx + 1] = (chi * step_val3 * inv0 + inv1) / delta;

    val0[0] = val1[0];
    val0[1] = val1[1];

    val1[0] = val2[0];
    val1[1] = val2[1];

    val2[0] = val3[0];
    val2[1] = val3[1];
  }
}

inline void lsms::rel_in(double energy, int k_qn, double *p_array, double *q_array,
                         const double *__restrict__ r_mesh,
                         double h,
                         const double *__restrict__ pot, std::size_t end,
                         double &p_last, double &q_last, std::size_t stop,
                         std::size_t &imax) {
  auto idx = end - 1;
  auto start_idx = idx;

  double val0[2], val1[2], val2[2];
  double delta, chi, inv0, inv1;

  double Mat11;
  double Mat21;
  double Mat12;
  double Mat22;

  constexpr auto tolerance = DBL_MIN;

  double omega = std::sqrt(-2 * energy - energy * energy / (C_SQ));

  
  for (idx = end - 1; idx >= stop; idx--) {
    p_array[idx] = std::exp(-omega * r_mesh[idx]) / std::sqrt(-energy / (energy + 2 * C_SQ));
    q_array[idx] = -std::exp(-omega * r_mesh[idx]);

    if (p_array[idx] > tolerance && std::abs(q_array[idx]) > tolerance) {
      start_idx = idx;
      break;
    }
  }

  imax = start_idx;

  p_array[start_idx - 1] =
      std::exp(-omega * r_mesh[start_idx - 1]) / std::sqrt(-energy / (energy + 2 * C_SQ));

  q_array[start_idx - 1] = -std::exp(-omega * r_mesh[start_idx - 1]);

  p_array[start_idx - 2] =
      std::exp(-omega * r_mesh[start_idx - 2]) / std::sqrt(-energy / (energy + 2 * C_SQ));

  q_array[start_idx - 2] = -std::exp(-omega * r_mesh[start_idx - 2]);

  idx = start_idx;
  rel_rad_func(r_mesh[idx], h * r_mesh[idx], energy, pot[idx], k_qn, p_array[idx], q_array[idx], val0);

  rel_rad_func(r_mesh[idx - 1], h * r_mesh[idx - 1], energy, pot[idx - 1], k_qn,
               p_array[idx - 1], q_array[idx - 1], val1);

  for (idx = start_idx; idx > stop + 3; --idx) {
    rel_rad_func(r_mesh[idx - 2], h * r_mesh[idx - 2], energy, pot[idx - 2], k_qn,
                 p_array[idx - 2], q_array[idx - 2], val2);

    inv0 = p_array[idx - 2] - 1.0 / 24 * (val0[0] - 5.0 * val1[0] + 19.0 * val2[0]);
    inv1 = q_array[idx - 2] - 1.0 / 24 * (val0[1] - 5.0 * val1[1] + 19.0 * val2[1]);
    chi = -9.0 / 24.0;

    delta = 1 + chi * chi * h * r_mesh[idx - 3] * h * r_mesh[idx - 3] *
        (-k_qn * k_qn / (r_mesh[idx - 3] * r_mesh[idx - 3]) -
            (-energy + pot[idx - 3]) * (2 * C + (energy - pot[idx - 3]) / C) / C);

    Mat11 = -k_qn * chi * h * r_mesh[idx - 3] / r_mesh[idx - 3] + 1;
    Mat22 = k_qn * chi * h * r_mesh[idx - 3] / r_mesh[idx - 3] + 1;
    Mat21 = chi * h * r_mesh[idx - 3] * (-energy + pot[idx - 3]) / C;
    Mat12 = chi * h * r_mesh[idx - 3] * (2 * C + (energy - pot[idx - 3]) / C);

    p_array[idx - 3] = (Mat11 * inv0 + Mat12 * inv1) / delta;
    q_array[idx - 3] = (Mat21 * inv0 + Mat22 * inv1) / delta;

    val0[0] = val1[0];
    val0[1] = val1[1];

    val1[0] = val2[0];
    val1[1] = val2[1];
  }

  idx = stop + 3;

  rel_rad_func(r_mesh[idx], h * r_mesh[idx], energy, pot[idx], k_qn, p_array[idx], q_array[idx], val0);

  rel_rad_func(r_mesh[idx - 1], h * r_mesh[idx - 1], energy, pot[idx - 1], k_qn,
               p_array[idx - 1], q_array[idx - 1], val1);

  rel_rad_func(r_mesh[idx - 2], h * r_mesh[idx - 2], energy, pot[idx - 2], k_qn,
               p_array[idx - 2], q_array[idx - 2], val2);

  inv0 = p_array[idx - 2] - 1.0 / 24.0 * (val0[0] - 5.0 * val1[0] + 19.0 * val2[0]);
  inv1 = q_array[idx - 2] - 1.0 / 24.0 * (val0[1] - 5.0 * val1[1] + 19.0 * val2[1]);

  chi = -9.0 / 24.0;

  delta = 1 + chi * chi * h * r_mesh[idx - 3] * h * r_mesh[idx - 3] *
      (-k_qn * k_qn / (r_mesh[idx - 3] * r_mesh[idx - 3]) -
          (-energy + pot[idx - 3]) * (2 * C + (energy - pot[idx - 3]) / C) / C);

  Mat11 = -k_qn * chi * h * r_mesh[idx - 3] / r_mesh[idx - 3] + 1;
  Mat22 = k_qn * chi * h * r_mesh[idx - 3] / r_mesh[idx - 3] + 1;

  Mat21 = chi * h * r_mesh[idx - 3] * (-energy + pot[idx - 3]) / C;
  Mat12 = chi * h * r_mesh[idx - 3] * (2 * C + (energy - pot[idx - 3]) / C);

  p_last = (Mat11 * inv0 + Mat12 * inv1) / delta;
  q_last = (Mat21 * inv0 + Mat22 * inv1) / delta;
}

inline void lsms::rel_out(double energy, int k_qn, double *p_array, double *q_array,
                          const double *__restrict__ r_mesh,
                          double h,
                          const double *__restrict__ pot, std::size_t end,
                          int sign_switch) {
  double yp[2], val0[2], val1[2], val2[2], val3[3];
  double delta, chi, inv0, inv1, Mat11, Mat21, Mat12, Mat22;

  double Z = -pot[0] * r_mesh[0];
  double beta = std::sqrt(k_qn * k_qn - (Z * Z) / (C * C));

  p_array[0] = std::pow(r_mesh[0], beta);
  q_array[0] = std::pow(r_mesh[0], beta) * C * (beta + k_qn) / Z;
  p_array[0] = p_array[0] * std::pow(-1, sign_switch);
  q_array[0] = q_array[0] * std::pow(-1, sign_switch);

  rel_rad_func(r_mesh[0], h * r_mesh[0], energy, pot[0], k_qn, p_array[0], q_array[0], val0);

  yp[0] = p_array[0] + val0[0];
  yp[1] = q_array[0] + val0[1];

  rel_rad_func(r_mesh[1], h * r_mesh[1], energy, pot[1], k_qn, yp[0], yp[1], val1);

  p_array[1] = p_array[0] + (val0[0] + val1[0]) / 2;
  q_array[1] = q_array[0] + (val0[1] + val1[1]) / 2;

  rel_rad_func(r_mesh[1], h * r_mesh[1], energy, pot[1], k_qn, p_array[1], q_array[1], val1);

  inv0 = p_array[1] + (-val0[0] + 8.0 * val1[0]) / 12.0;
  inv1 = q_array[1] + (-val0[1] + 8.0 * val1[1]) / 12.0;

  chi = 5.0 / 12.0;

  delta = 1 + chi * chi * h * r_mesh[2] * h * r_mesh[2] *
      (-k_qn * k_qn / (r_mesh[2] * r_mesh[2]) -
          (-energy + pot[2]) * (2 * C + (energy - pot[2]) / C) / C);

  Mat11 = -k_qn * chi * h * r_mesh[2] / r_mesh[2] + 1;
  Mat22 = k_qn * chi * h * r_mesh[2] / r_mesh[2] + 1;
  Mat21 = chi * h * r_mesh[2] * (-energy + pot[2]) / C;
  Mat12 = chi * h * r_mesh[2] * (2 * C + (energy - pot[2]) / C);

  p_array[2] = (Mat11 * inv0 + Mat12 * inv1) / delta;
  q_array[2] = (Mat21 * inv0 + Mat22 * inv1) / delta;

  rel_rad_func(r_mesh[2], h * r_mesh[2], energy, pot[2], k_qn, p_array[2], q_array[2], val2);

  inv0 = p_array[2] + (val0[0] - 5.0 * val1[0] + 19.0 * val2[0]) / 24.0;
  inv1 = q_array[2] + (val0[1] - 5.0 * val1[1] + 19.0 * val2[1]) / 24.0;

  chi = 9.0 / 24.0;
  delta = 1 + chi * chi * h * r_mesh[3] * h * r_mesh[3] *
      (-k_qn * k_qn / (r_mesh[3] * r_mesh[3]) -
          (-energy + pot[3]) * (2 * C + (energy - pot[3]) / C) / C);

  Mat11 = -k_qn * chi * h * r_mesh[3] / r_mesh[3] + 1;
  Mat22 = k_qn * chi * h * r_mesh[3] / r_mesh[3] + 1;
  Mat21 = chi * h * r_mesh[3] * (-energy + pot[3]) / C;
  Mat12 = chi * h * r_mesh[3] * (2 * C + (energy - pot[3]) / C);

  p_array[3] = (Mat11 * inv0 + Mat12 * inv1) / delta;
  q_array[3] = (Mat21 * inv0 + Mat22 * inv1) / delta;

  for (int idx = 3; idx < end - 1; ++idx) {
    rel_rad_func(r_mesh[idx], h * r_mesh[idx], energy, pot[idx], k_qn, p_array[idx], q_array[idx],
                 val3);

    inv0 = p_array[idx] +
        1.0 / 720 *
            (-19.0 * val0[0] + 106.0 * val1[0] - 264.0 * val2[0] + 646.0 * val3[0]);

    inv1 = q_array[idx] +
        1.0 / 720 *
            (-19.0 * val0[1] + 106.0 * val1[1] - 264.0 * val2[1] + 646.0 * val3[1]);

    chi = 251.0 / 720.0;

    delta = 1 + chi * chi * h * r_mesh[idx + 1] * h * r_mesh[idx + 1] *
        (-k_qn * k_qn / (r_mesh[idx + 1] * r_mesh[idx + 1]) -
            (-energy + pot[idx + 1]) * (2 * C + (energy - pot[idx + 1]) / C) / C);

    Mat11 = -k_qn * chi * h * r_mesh[idx + 1] / r_mesh[idx + 1] + 1;
    Mat22 = k_qn * chi * h * r_mesh[idx + 1] / r_mesh[idx + 1] + 1;
    Mat21 = chi * h * r_mesh[idx + 1] * (-energy + pot[idx + 1]) / C;
    Mat12 = chi * h * r_mesh[idx + 1] * (2 * C + (energy - pot[idx + 1]) / C);

    p_array[idx + 1] = (Mat11 * inv0 + Mat12 * inv1) / delta;

    q_array[idx + 1] = (Mat21 * inv0 + Mat22 * inv1) / delta;

    val0[0] = val1[0];
    val0[1] = val1[1];

    val1[0] = val2[0];
    val1[1] = val2[1];

    val2[0] = val3[0];
    val2[1] = val3[1];
  }
}

inline double lsms::rel_integration(const double *__restrict__ p_array,
                                    const double *__restrict__ q_array,
                                    const double *__restrict__ r_mesh,
                                    double h, std::size_t end) {
  double cummulative_val = 0.0;

  int idx;

  for (idx = 0; idx < end - 2; idx += 2) {
    cummulative_val +=
        1.0 / 3.0 *
            ((p_array[idx] * p_array[idx] + q_array[idx] * q_array[idx]) * h * r_mesh[idx] +
                4 * (p_array[idx + 1] * p_array[idx + 1] + q_array[idx + 1] * q_array[idx + 1]) * h * r_mesh[idx + 1] +
                (p_array[idx + 2] * p_array[idx + 2] + q_array[idx + 2] * q_array[idx + 2]) * h * r_mesh[idx + 2]);
  }

  if ((end % 2) == 0) {
    cummulative_val += 0.5 * ((p_array[idx] * p_array[idx] + q_array[idx] * q_array[idx]) * h * r_mesh[idx] +
        (p_array[idx + 1] * p_array[idx + 1] + q_array[idx + 1] * q_array[idx + 1]) *
            h * r_mesh[idx + 1]);
  }

  return cummulative_val;
}

inline double lsms::nrel_integration(const double *__restrict__ P,
                                     const double *__restrict__ r_mesh,
                                     double h,
                                     std::size_t end) {
  double cummulative_val = 0.0;

  int idx;

  for (idx = 0; idx < end - 2; idx += 2) {
    cummulative_val += 1.0 / 3.0 *
        ((P[idx] * P[idx]) * h * r_mesh[idx] +
            4 * (P[idx + 1] * P[idx + 1]) * h * r_mesh[idx + 1] +
            (P[idx + 2] * P[idx + 2]) * h * r_mesh[idx + 2]);
  }

  if ((end % 2) == 0) {
    cummulative_val += 0.5 * ((P[idx] * P[idx]) * h * r_mesh[idx] +
        (P[idx + 1] * P[idx + 1]) * h * r_mesh[idx + 1]);
  }

  return cummulative_val;
}

double lsms::rel_eigenenergies(
    double Z, int n_qn, int l_qn, int k_qn, double *dens, double *p_array, double *q_array,
    const double *__restrict__ r_mesh, double h,
    const double *__restrict__ pot, std::size_t end, double energy_rel_tol, int max_iter,
    double e_init, int &converged, double &delta_energy) {

  double emax_init = 1.0e-6;
  double emin_init = -10000000;

  if (l_qn != 0) {
    emin_init = pot[0] + (l_qn * (l_qn + 1)) / (2.0 * r_mesh[0] * r_mesh[0]);

    for (int ir = 1; ir < end; ir++) {
      emin_init =
          std::min((pot[ir] + l_qn * (l_qn + 1) / (2.0 * r_mesh[ir] * r_mesh[ir])), emin_init);
    }

    emin_init = emin_init * 1.1;

  } else {
    emin_init = rel_energy_start(n_qn, k_qn, Z) * 1.6;
  }

  if (e_init == 0.0) {
    e_init = -(Z * Z) / ((2.0 * n_qn * n_qn));
  } else if (e_init < emin_init || emax_init > 0.0) {
    e_init = emin_init / 4.0;
  } else {
    e_init = e_init;
  }

  double E = e_init;

  const int number_of_nodes = n_qn - 1 - l_qn;
  int counted_counted_nr_nodes;

  auto i = max_iter;
  auto j = end;
  auto ctp = end;

  double dE = 0.0;
  double P_minus, P_plus, Q_minus, Q_plus;

  double factor;
  double norm_val;
  int sign_switch;

  int iter_max = max_iter;
  auto last_point = end;
  double etol = energy_rel_tol;

  double e_max = emax_init;
  double e_min = emin_init;

  bool last_bisect = true;
  bool too_many_switches = true;
  bool false_nr_nodes = true;

  for (i = 0; i < iter_max; i++) {
    
    ctp = calculate_ctp(E, l_qn, r_mesh, pot, end);

    if (ctp > end - 10) {
      ctp = end - 1;
    }

    if (E >= 0) {
      ctp = end - 1;
    }

    sign_switch = number_of_nodes;
    rel_out(E, k_qn, p_array, q_array, r_mesh, h, pot, ctp + 1, sign_switch);

    auto imax = last_valid_point(p_array, end);

    if (imax >= ctp) {
      imax = ctp;
    }

    counted_counted_nr_nodes = count_sign_sign_switch(p_array, 0, imax + 1);

    if (number_of_nodes != counted_counted_nr_nodes || ctp == end - 1 ||
        imax < ctp) {
      if (counted_counted_nr_nodes > number_of_nodes) {
        too_many_switches = true;
      } else {
        too_many_switches = false;
      }

      if (too_many_switches) {
        e_max = E;
      } else {
        e_min = E;
      }

      E = (e_min + e_max) / 2;

      last_bisect = true;
      false_nr_nodes = true;
      continue;
    } else {
      false_nr_nodes = false;
    }

    rel_in(E, k_qn, p_array, q_array, r_mesh, h, pot, end, P_plus, Q_plus, ctp,
           last_point);

    P_minus = p_array[ctp];
    Q_minus = q_array[ctp];

    factor = P_minus / P_plus;

    for (j = ctp + 1; j < end; j++) {
      p_array[j] = p_array[j] * factor;
      q_array[j] = q_array[j] * factor;
    }

    Q_plus = Q_plus * factor;
    P_plus = P_plus * factor;

    norm_val = rel_integration(p_array, q_array, r_mesh, h, last_point);
    dE = C * (Q_minus - Q_plus) * P_minus / (norm_val);

    if (std::fabs(dE / E) < etol) {
      double sq_norm = std::sqrt(norm_val);

      for (j = 0; j < end; j++) {
        q_array[j] = q_array[j] / sq_norm;
        p_array[j] = p_array[j] / sq_norm;
        dens[j] = p_array[j] * p_array[j] + q_array[j] * q_array[j];
      }

      converged = i;
      delta_energy = std::fabs(dE / E);
      return E;
    }

    too_many_switches = dE < 0;

    if (too_many_switches) {
      e_max = E;
    } else {
      e_min = E;
    }

    if (E + dE > e_max || E + dE < e_min) {
      E = (e_min + e_max) / 2.0;
      last_bisect = true;
    } else {
      E = E + dE;
      last_bisect = false;
    }
  }

  converged = ERROR_NOT_REACHED_ACCURACY;
  delta_energy = std::fabs(dE / E);

  if (last_bisect) {
    converged = ERROR_STILL_BISECTION;
  }

  if (false_nr_nodes) {
    converged = ERROR_WRONG_NUMBER_NODES;
  }

  norm_val = rel_integration(p_array, q_array, r_mesh, h, end);

  double sq_norm = std::sqrt(norm_val);

  for (j = 0; j < end; j++) {
    q_array[j] = q_array[j] / sq_norm;
    p_array[j] = p_array[j] / sq_norm;
    dens[j] = p_array[j] * p_array[j] + q_array[j] * q_array[j];
  }

  return E;
}

double lsms::nonrel_eigenenergies(
    double Z, int n_qn, int l_qn, int k_qn, double *dens, double *p_array, double *q_array,
    const double *__restrict__ r_mesh, double h,
    const double *__restrict__ pot, std::size_t end, double energy_rel_tol, int max_iter,
    double e_init, int &converged, double &delta_energy) {
  double emax_init;
  double emin_init;

  
  emax_init = 1.0e-6;

  if (l_qn != 0) {
    emin_init = pot[0] + (l_qn * (l_qn + 1)) / (2.0 * r_mesh[0] * r_mesh[0]);

    for (int ir = 1; ir < end; ir++) {
      emin_init =
          std::min((pot[ir] + l_qn * (l_qn + 1) / (2.0 * r_mesh[ir] * r_mesh[ir])), emin_init);
    }

    emin_init = emin_init * 1.1;

  } else {
    emin_init = nrel_energy_start(n_qn, Z) * 1.6;
  }

  
  if (e_init == 0.0) {
    e_init = -(Z * Z) / ((3.0 * n_qn * n_qn));
  } else if (e_init < emin_init || emax_init > 0.0) {
    e_init = emin_init / 3.0;
  } else {
    e_init = e_init;
  }

  double E = e_init;

  
  const int number_of_nodes = n_qn - 1 - l_qn;

  int counted_counted_nr_nodes;
  auto i = max_iter;
  auto j = end;
  auto ctp = end;
  auto last_point = end;

  double dE = 0.0;
  double P_minus, P_plus, Q_minus, Q_plus;
  double factor, norm_val;

  int sign_switch;
  double etol = energy_rel_tol;

  double e_max = emax_init;
  double e_min = emin_init;

  bool last_bisect = true;
  bool too_many_switches = true;
  bool false_nr_nodes = true;

  for (i = 0; i < max_iter; i++) {
    
    ctp = calculate_ctp(E, l_qn, r_mesh, pot, end);

    if (ctp > end - 10) {
      ctp = end - 1;
    }

    if (E >= 0) {
      ctp = end - 1;
    }

    sign_switch = number_of_nodes;
    nrel_out(Z, E, l_qn, p_array, q_array, r_mesh, h, pot, ctp + 1, sign_switch);

    auto imax = last_valid_point(p_array, end);

    if (imax >= ctp) {
      imax = ctp;
    }

    counted_counted_nr_nodes = count_sign_sign_switch(p_array, 0, imax + 1);

    if (number_of_nodes != counted_counted_nr_nodes || ctp == end - 1 ||
        imax < ctp) {
      if (counted_counted_nr_nodes > number_of_nodes) {
        too_many_switches = true;
      } else {
        too_many_switches = false;
      }

      if (too_many_switches) {
        e_max = E;
      } else {
        e_min = E;
      }

      E = (e_min + e_max) / 2;
      last_bisect = true;
      false_nr_nodes = true;
      continue;
    } else {
      false_nr_nodes = false;
    }

    nrel_in(E, l_qn, p_array, q_array, r_mesh, h, pot, end, P_plus, Q_plus, ctp, last_point);

    P_minus = p_array[ctp];
    Q_minus = q_array[ctp];
    factor = P_minus / P_plus;

    for (j = ctp + 1; j < end; j++) {
      p_array[j] = p_array[j] * factor;
      q_array[j] = q_array[j] * factor;
    }

    Q_plus = Q_plus * factor;
    P_plus = P_plus * factor;

    norm_val = nrel_integration(p_array, r_mesh, h, last_point);

    dE = (Q_minus - Q_plus) * P_minus / (2 * norm_val);

    if (std::fabs(dE / E) < etol) {
      double sq_norm = std::sqrt(norm_val);

      for (j = 0; j < end; j++) {
        q_array[j] = q_array[j] / sq_norm;
        p_array[j] = p_array[j] / sq_norm;
        dens[j] = p_array[j] * p_array[j];
      }

      converged = i;
      delta_energy = std::fabs(dE / E);
      return E;
    }

    
    too_many_switches = dE < 0;

    if (too_many_switches) {
      e_max = E;
    } else {
      e_min = E;
    }

    if (E + dE > e_max || E + dE < e_min) {
      E = (e_min + e_max) / 2.0;
      last_bisect = true;
    } else {
      E = E + dE;
      last_bisect = false;
    }
  }

  converged = ERROR_NOT_REACHED_ACCURACY;
  delta_energy = std::fabs(dE / E);

  if (last_bisect) {
    converged = ERROR_STILL_BISECTION;
  }

  if (false_nr_nodes) {
    converged = ERROR_WRONG_NUMBER_NODES;
  }

  norm_val = nrel_integration(p_array, r_mesh, h, end);
  auto sq_norm = std::sqrt(norm_val);

  for (j = 0; j < end; j++) {
    q_array[j] = q_array[j] / sq_norm;
    p_array[j] = p_array[j] / sq_norm;
    dens[j] = p_array[j] * p_array[j];
  }

  return E;
}

double lsms::rel_energy_start(int n, int k_qn, double Z) {
  auto beta = std::sqrt(k_qn * k_qn - Z * Z / (C * C));

  return C * C /
      std::sqrt(1 + (Z * Z / (C * C)) /
          ((n - abs(k_qn) + beta) * (n - abs(k_qn) + beta))) -
      C * C;
}

double lsms::nrel_energy_start(int n, double Z) { return -Z * Z / (2.0 * n * n); }
