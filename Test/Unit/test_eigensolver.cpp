//
// Created by F.Moitzi on 15.05.2022.
//

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>

#include "Core/States.hpp"
#include "Core/radialSolver.hpp"
#include "Misc/diff.hpp"
#include "PhysicalConstants.hpp"
#include "accel_common.hpp"

namespace EigenSolverTests {

template <typename T>
T relative_diff(T ref, T val) {
  return std::fabs(ref - val) / std::fabs(ref);
}

extern "C" {
void deepst_(int *nqn, int *lqn, int *kqn, double *en, double *rv, double *r,
             double *rf, double *h, double *z, double *c, int *nitmax,
             double *tol, int *nws, int *nlast, int *iter, int *iprpts,
             int *ipdeq);
}

TEST(EigenSolverTests, NonRelTest) {
  int N = 2000;
  std::vector<double> dens(N);
  std::vector<double> rf(N);
  std::vector<double> R(N);
  std::vector<double> Rp(N);
  std::vector<double> P(N);
  std::vector<double> Q(N);
  std::vector<double> V(N);
  std::vector<double> vr(N);

  double r0 = 0.001;
  double h = 0.03;
  double Z = 92;
  r0 /= Z;

  // New solver
  double energy_rel_tol = 1.0e-15;
  int max_iter = 150;
  double e_init = 0.0;
  int converged;
  double delta_energy;

  // Legacy solver
  double c = cphot * 100000000;
  int nitmax = 300;
  int ipdeq = 5;
  int iter;
  double tol = 1.0e-10;

  for (auto ir{0}; ir < N; ir++) {
    R[ir] = r0 * std::exp(h * ir);
    Rp[ir] = R[ir] * h;
    V[ir] = -Z / R[ir];
    vr[ir] = V[ir] * 2.0 * R[ir];
  }

  double Rws = 3.5;
  int jws = std::floor(std::log(Rws / r0) / h);
  Rws = R[jws - 1];
  std::cout << Rws << std::endl;

  std::vector<int> n_states;
  std::vector<int> l_states;
  std::vector<int> spin_states;
  std::vector<int> kappa_states;
  std::vector<double> occ_states;

  lsms::States::relativistic_atomic_states(Z, n_states, l_states, spin_states,
                                           kappa_states, occ_states);

  for (int index = 0; index < n_states.size(); index++) {
    int n = n_states[index];
    int l = l_states[index];
    int spin = spin_states[index];
    int kappa = kappa_states[index];

    std::cout << " ----------- " << std::endl;
    std::cout << n << " " << l << " " << kappa << std::endl;

    double energy = lsms::nonrel_eigenenergies(
        Z, n, l, 0, dens.data(), P.data(), Q.data(), R.data(), h, V.data(), N,
        energy_rel_tol, max_iter, e_init, converged, delta_energy);

    ASSERT_TRUE(converged >= 0);

    double energy_ref = -Z * Z / (2.0 * n * n);

    // Output
    // std::cout << energy << std::endl;
    std::cout << relative_diff(energy, energy_ref) << std::endl;
    ASSERT_LE(relative_diff(energy, energy_ref), 10e-6);

    /***
     * Reference
     */

    double energy_leg = energy_ref / 4.0;

    std::fill(rf.begin(), rf.end(), 0.0);

    deepst_(&n, &l, &kappa, &energy_leg, vr.data(), R.data(), rf.data(), &h, &Z,
            &c, &nitmax, &tol, &jws, &N, &iter, &N, &ipdeq);

    // std::cout << energy_leg / 2.0 << std::endl;
    std::cout << relative_diff(energy_leg / 2.0, energy_ref) << std::endl;
  }
}

TEST(EigenSolverTests, RelTest) {
  int N = 2000;
  std::vector<double> dens(N);
  std::vector<double> rf(N);
  std::vector<double> R(N);
  std::vector<double> Rp(N);
  std::vector<double> P(N);
  std::vector<double> Q(N);
  std::vector<double> V(N);
  std::vector<double> vr(N);

  double r0 = 0.001;
  double h = 0.03;
  double Z = 92;
  r0 /= Z;

  // New solver
  double energy_rel_tol = 1.0e-10;
  int max_iter = 300;
  double e_init = 0.0;
  int converged;
  double delta_energy;

  // Legacy solver
  double c = cphot;
  int nitmax = 300;
  int ipdeq = 5;
  int iter;
  double tol = 1.0e-10;

  for (auto ir{0}; ir < N; ir++) {
    R[ir] = r0 * std::exp(h * ir);
    Rp[ir] = R[ir] * h;
    V[ir] = -Z / R[ir];
    vr[ir] = V[ir] * 2.0 * R[ir];
  }

  double Rws = 3.5;
  int jws = std::floor(std::log(Rws / r0) / h);
  Rws = R[jws - 1];
  std::cout << Rws << std::endl;

  std::vector<double> refs = {
      -4861.19790417167860141490, -1257.39585206415063112217,
      -1257.39585206368451508752, -1089.61141621799083623046,
      -539.09332897057436184696,  -539.09332896996170347848,
      -489.03708486665664167958,  -489.03708486600510241260,
      -476.26159429296740199788,  -295.25783850197689162087,
      -295.25783850149235831850,  -274.40775736073044299701,
      -274.40775735994373007998,  -268.96587718349456963551,
      -268.96587718285354640102,  -266.38944691927628127814,
      -185.48518877602776910862,  -185.48518877565157936260,
      -174.94461273937429268699,  -174.94461273871903017607,
      -172.15525190859702320267,  -172.15525190779592890067,
      -170.82893682911284827242,  -127.09363715385616444564,
      -127.09363715353399015839,  -121.05753750855492967275,
      -121.05753750798186274551,  -119.44527171386766895012,
      -92.44078653463942885082};

  std::vector<int> n_states;
  std::vector<int> l_states;
  std::vector<int> spin_states;
  std::vector<int> kappa_states;
  std::vector<double> occ_states;

  lsms::States::relativistic_atomic_states(Z, n_states, l_states, spin_states,
                                           kappa_states, occ_states);

  for (int index = 0; index < n_states.size(); index++) {
    int n = n_states[index];
    int l = l_states[index];
    int spin = spin_states[index];
    int kappa = kappa_states[index];

    std::cout << " ----------- " << std::endl;
    std::cout << n << " " << l << " " << kappa << std::endl;

    double energy = lsms::rel_eigenenergies(
        Z, n, l, kappa, dens.data(), P.data(), Q.data(), R.data(), h, V.data(),
        N, energy_rel_tol, max_iter, e_init, converged, delta_energy);

    ASSERT_TRUE(converged >= 0);

    double energy_ref = -Z * Z / (2.0 * n * n);

    // Output
    std::cout << energy << std::endl;
    std::cout << relative_diff(energy, refs[index]) << std::endl;
    ASSERT_LE(relative_diff(energy, refs[index]), 10e-6);

    /***
     * Reference
     */

    double energy_leg = energy_ref / 4.0;

    std::fill(rf.begin(), rf.end(), 0.0);

    deepst_(&n, &l, &kappa, &energy_leg, vr.data(), R.data(), rf.data(), &h, &Z,
            &c, &nitmax, &tol, &jws, &N, &iter, &N, &ipdeq);

    std::cout << energy_leg / 2.0 << std::endl;
    std::cout << relative_diff(energy_leg / 2.0, refs[index]) << std::endl;
  }
}

}  // namespace EigenSolverTests