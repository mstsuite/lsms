//
// Created by F.Moitzi on 15.05.2022.
//

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <numeric>

#include "Core/States.hpp"
#include "Core/atomic_dft.hpp"
#include "Mixer.hpp"
#include "MixingVector.hpp"
#include "XCLDA.hpp"
#include "accel_common.hpp"
#include "fmt/core.h"
#include "fmt/ranges.h"
#include "poisson.hpp"
#include "xc_lda_potential.hpp"

namespace atomic_dft {

class AtomicDFTTestFixture
    : public ::testing::TestWithParam<std::tuple<int, Real>> {};

TEST(AtomicDFTTest, Rmesh) {
  Real r_min = 0.001;  // starting of rmesh
  Real r_max = 50.0;   // ending of rmesh
  Real h_step = 0.02;  // default h step in atom

  Real ws = 2.923;

  int max_point = std::floor(std::log(r_max / r_min) / h_step) + 1;

  int jmt = std::floor(std::log(ws / r_min) / h_step) + 1;

  Real rmin_new = ws / std::exp(h_step * (jmt - 1));

  std::cout << rmin_new << std::endl;

  std::cout << std::log(rmin_new) << std::endl;

  ASSERT_NEAR(rmin_new * std::exp(h_step * (jmt - 1)), ws, 1.0e-8);
}

TEST_P(AtomicDFTTestFixture, RLDA) {
  int Z = std::get<0>(GetParam());
  Real tot_energy_ref = std::get<1>(GetParam());

  lsms::AtomicDFT dft;

  int max_iter = 400;
  int max_eig_iter = 100;
  Real e_eig_tol = 1.0e-12;
  Real e_tol = 1.0e-9;

  int N = 1300;
  Real r0 = 0.0001;
  Real h = 0.01;

  std::vector<Real> r_mesh(N);

  for (auto ir{0}; ir < N; ir++) {
    r_mesh[ir] = r0 * std::exp(h * ir);
  }

  std::vector<Real> tot_potential(N);
  Matrix<Real> density(N, 1);

  auto [ks_energies, tot_energy] =
      dft.solve(Z, r_mesh, h, N, density, tot_potential);
  EXPECT_LE(
      std::fabs(tot_energy - 2.0 * tot_energy_ref) / std::fabs(tot_energy),
      5.0e-6);
}

TEST_P(AtomicDFTTestFixture, LibxcVWM) {
  int Z = std::get<0>(GetParam());
  Real tot_energy_ref = std::get<1>(GetParam());

  std::vector<int> functionals{1, 532, 7};

  lsms::AtomicDFT dft(functionals);

  int max_iter = 400;
  int max_eig_iter = 100;
  Real e_eig_tol = 1.0e-12;
  Real e_tol = 1.0e-9;

  int N = 1300;
  Real r0 = 0.0001;
  Real h = 0.01;

  std::vector<Real> r_mesh(N);

  for (auto ir{0}; ir < N; ir++) {
    r_mesh[ir] = r0 * std::exp(h * ir);
  }

  std::vector<Real> tot_potential(N);
  Matrix<Real> density(N, 1);

  auto [ks_energies, tot_energy] =
      dft.solve(Z, r_mesh, h, N, density, tot_potential);
  EXPECT_LE(
      std::fabs(tot_energy - 2.0 * tot_energy_ref) / std::fabs(tot_energy),
      5.0e-6);
}

TEST_P(AtomicDFTTestFixture, LibxcGGAPW91) {
  int Z = std::get<0>(GetParam());

  std::vector<int> functionals{1, 109, 134};

  int max_iter = 400;
  int max_eig_iter = 100;
  Real e_eig_tol = 1.0e-13;
  Real e_tol = 1.0e-10;

  int N = 1300;
  Real r0 = 0.0001;
  Real h = 0.01;

  lsms::AtomicDFT dft(functionals, max_iter, max_eig_iter, e_tol, e_eig_tol);

  std::vector<Real> r_mesh(N);

  for (auto ir{0}; ir < N; ir++) {
    r_mesh[ir] = r0 * std::exp(h * ir);
  }

  std::vector<Real> tot_potential(N);
  Matrix<Real> density(N, 1);

  auto [ks_energies, tot_energy] =
      dft.solve(Z, r_mesh, h, N, density, tot_potential);
}

TEST_P(AtomicDFTTestFixture, TotalEnergies) {
  int Z = std::get<0>(GetParam());
  Real tot_energy_ref = std::get<1>(GetParam());

  int max_iter = 400;
  int max_eig_iter = 100;
  Real e_eig_tol = 1.0e-12;
  Real e_tol = 1.0e-9;

  int N = 1300;
  Real r0 = 0.0001;
  Real h = 0.01;

  std::vector<Real> potential(N);
  std::vector<Real> tot_potential(N);
  std::vector<Real> potential_prev(N);

  std::vector<Real> vhartree_potential(N);
  std::vector<Real> vhartree_potential_ref(N);
  std::vector<Real> dvhartree_potential(N);
  std::vector<Real> r_mesh(N);

  Matrix<Real> xc_potential(N, 1);
  Matrix<Real> e_xc(N, 1);
  Matrix<Real> density(N, 1);
  Matrix<Real> density_prev(N, 1);

  for (auto ir{0}; ir < N; ir++) {
    r_mesh[ir] = r0 * std::exp(h * ir);
  }

  std::vector<int> n;
  std::vector<int> l;
  std::vector<int> spin;
  std::vector<int> kappa;
  std::vector<Real> occupation;

  lsms::DefaultMixerVector mixVector(N);

  lsms::XCLDA xc_lda(1, {0, 0, 0});

  std::unique_ptr<lsms::AbstractMixer> mixer =
      std::make_unique<lsms::BroydenMixer>(mixVector.total_size, 0.05, 10, 25,
                                           0.01);

  auto start = std::chrono::high_resolution_clock::now();

  lsms::States::relativistic_atomic_states(Z, n, l, spin, kappa, occupation);

  std::vector<Real> e_ks_eigenenergies(n.size(), 0.0);

  fmt::print("n:    {:4d}\n", fmt::join(n, " "));
  fmt::print("l:    {:4d}\n", fmt::join(l, " "));
  fmt::print("spin: {:4d}\n", fmt::join(spin, " "));
  fmt::print("occ:  {:4.2f}\n", fmt::join(occupation, " "));

  Real total =
      std::accumulate(occupation.begin(), occupation.end(), 0.0, std::plus<>());

  ASSERT_NEAR(total, Z, 1.0e-12);

  // Generate starting potential
  lsms::generate_starting_potential(tot_potential, r_mesh, Z, N);

  for (auto ir = 0; ir < N; ir++) {
    mixVector.vec_new[ir] = tot_potential[ir] + Z / r_mesh[ir];
    mixVector.vec_old[ir] = tot_potential[ir] + Z / r_mesh[ir];
  }

  /*
   * SCF iteration
   */

  Real energy;
  Real tot_energy_prev = 0.0;
  Real tot_energy = 0.0;
  int iter = 0;

  for (iter = 0; iter < max_iter; iter++) {
    for (auto ir = 0; ir < N; ir++) {
      tot_potential[ir] = -Z / r_mesh[ir] + mixVector.vec_new[ir];
    }

    // Generate density
    energy = lsms::generate_density(r_mesh, h, N, tot_potential, Z, n, l, spin,
                                    kappa, occupation, e_ks_eigenenergies,
                                    density, max_eig_iter, e_eig_tol);

    // XC potential
    xc_lda.evaluate(r_mesh, h, density, N, e_xc, xc_potential);

    // Hartree potential
    lsms::radial_poisson(vhartree_potential, dvhartree_potential, r_mesh, h,
                         density, N);

    // Total energy
    tot_energy = lsms::total_energy(r_mesh, density, vhartree_potential, e_xc,
                                    tot_potential, energy, Z, N);

    // Mix density
    for (auto ir = 0; ir < N; ir++) {
      mixVector.vec_new[ir] = xc_potential[ir] * 0.5 + vhartree_potential[ir];
    }

    Real vrms = 0.0;
    Real qrms = 0.0;
    for (auto ir = 0; ir < N; ir++) {
      vrms += (mixVector.vec_new[ir] - mixVector.vec_old[ir]) *
              (mixVector.vec_new[ir] - mixVector.vec_old[ir]);
      qrms +=
          (density[ir] - density_prev[ir]) * (density[ir] - density_prev[ir]);
    }

    mixVector.rms = std::sqrt(vrms);
    qrms = std::sqrt(qrms);

    mixer->mix(mixVector);

    for (auto ir = 0; ir < N; ir++) {
      mixVector.vec_old[ir] = mixVector.vec_new[ir];
      density_prev[ir] = density[ir];
    }

    if (std::abs(tot_energy - tot_energy_prev) < e_tol) {
      fmt::print(
          "{:4d} Band energy: {:17.10f} Energy: {:17.10f} VRMS: {:10.4e} QRMS: "
          "{:10.4e}\n",
          iter, energy, tot_energy, vrms, qrms);
      break;
    }

    tot_energy_prev = tot_energy;
  }

  auto end = std::chrono::high_resolution_clock::now();

  auto duration =
      std::chrono::duration_cast<std::chrono::microseconds>(end - start);

  fmt::print("{}\n", duration.count());

  EXPECT_LE(std::fabs(tot_energy_prev - 2.0 * tot_energy_ref) /
                std::fabs(tot_energy_prev),
            5.0e-6);
}

INSTANTIATE_TEST_SUITE_P(
    AtomicDFTTests, AtomicDFTTestFixture,
    ::testing::Values(std::make_tuple(1, -0.44566816250964969015),
                      std::make_tuple(3, -7.33523068192602778481),
                      std::make_tuple(11, -161.60169069548197740005),
                      std::make_tuple(13, -241.66932269398179755626),
                      std::make_tuple(91, -27223.30304328585771145299)));

}  // namespace atomic_dft