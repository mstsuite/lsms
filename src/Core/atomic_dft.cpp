//
// Created by F.Moitzi on 10.01.2023.
//

#include "atomic_dft.hpp"

#include <cmath>
#include <iostream>
#include <numeric>

#include "fmt/core.h"
#include "fmt/ranges.h"

#include "Core/States.hpp"
#include "Core/atomic_dft.hpp"

#include "xc_lda_potential.hpp"
#include "poisson.hpp"
#include "Mixer.hpp"
#include "MixingVector.hpp"
#include "radialSolver.hpp"
#include "integrator.hpp"
#include "XCLibxc.hpp"

void lsms::generate_starting_potential(std::vector<Real> &potential,
                                       const std::vector<Real> &r_mesh,
                                       int core_charge,
                                       std::size_t end
) {

  int Z = core_charge;

  std::vector<Real> Z_eff(end);
  std::vector<Real> x(end);

  for (int i = 0; i < end; i++) {
    x[i] = r_mesh[i] * std::pow((128 * Z / (9 * M_PI * M_PI)), 1.0 / 3.0);
  }

  Real alpha = 0.7280642371;
  Real beta = -0.5430794693;
  Real gamma = 0.3612163121;

  for (int i = 0; i < end; i++) {
    Z_eff[i] = Z *
        (1 + alpha * std::sqrt(x[i]) + beta * x[i] * exp(-gamma * std::sqrt(x[i]))) *
        (1 + alpha * std::sqrt(x[i]) + beta * x[i] * exp(-gamma * std::sqrt(x[i]))) *
        exp(-2 * alpha * std::sqrt(x[i]));

    if (Z_eff[i] < 1) {
      Z_eff[i] = 1;
    }

    potential[i] = -Z_eff[i] / r_mesh[i];
  }

}

Real lsms::generate_density(
    const std::vector<Real> &r_mesh,
    Real h,
    std::size_t end,
    const std::vector<Real> &potential,
    int core_charge,
    const std::vector<int> &n,
    const std::vector<int> &l,
    const std::vector<int> &spin,
    const std::vector<int> &kappa,
    const std::vector<Real> &occupation,
    std::vector<Real> &e_eig,
    Matrix<Real> &density,
    int max_eig_iter,
    Real eig_tol
) {

  std::vector<Real> dens(end);
  std::vector<Real> P(end);
  std::vector<Real> Q(end);

  auto max_orbs = n.size();

  Real energy = 0.0;

  density = 0.0;

  for (auto i_orbs = 0; i_orbs < max_orbs; i_orbs++) {

    Real e_init = e_eig[i_orbs];

    if (e_init == 0.0) {
      e_init = rel_energy_start(n[i_orbs], kappa[i_orbs], core_charge);
    }

    int converged = 0;
    Real delta_energy = 0.0;

    auto eig_energy = lsms::rel_eigenenergies(
        core_charge,
        n[i_orbs],
        l[i_orbs],
        kappa[i_orbs],
        dens.data(),
        P.data(),
        Q.data(),
        r_mesh.data(),
        h,
        potential.data(),
        end,
        eig_tol,
        max_eig_iter,
        e_init,
        converged,
        delta_energy) * occupation[i_orbs];

    e_eig[i_orbs] = eig_energy;

    energy += eig_energy;

    for (auto ir = 0; ir < end; ir++) {
      density(ir, 0) += dens[ir] * occupation[i_orbs];
    }

  }

  return energy;

}

Real lsms::total_energy(const std::vector<Real> &r_mesh,
                        Matrix<Real> &rho,
                        std::vector<Real> &v_hartree,
                        Matrix<Real> &e_xc,
                        std::vector<Real> &v_pot,
                        Real E_band,
                        Real Z,
                        int N
) {

  std::vector<Real> integrand(N);

  for (int ir = 0; ir < N; ir++) {
    integrand[ir] = v_pot[ir] * rho(ir, 0);
  }

  Real T_s = E_band - lsms::radialIntegral(integrand, r_mesh, N);

  for (int ir = 0; ir < N; ir++) {
    integrand[ir] = -0.5 * v_hartree[ir] * rho(ir, 0);
  }

  Real E_ee = -lsms::radialIntegral(integrand, r_mesh, N);

  for (int ir = 0; ir < N; ir++) {
    integrand[ir] = -Z / r_mesh[ir] * rho(ir, 0);
  }

  Real E_en = lsms::radialIntegral(integrand, r_mesh, N);
  Real E_c = E_ee + E_en;

  for (int ir = 0; ir < N; ir++) {
    integrand[ir] = e_xc(ir, 0) * rho(ir, 0);
  }

  Real EE_xc = lsms::radialIntegral(integrand, r_mesh, N);

  return 2.0 * T_s + 2.0 * E_c + EE_xc;
}

std::tuple<std::vector<double>, double> lsms::AtomicDFT::solve(int Z,
                                                               const std::vector<Real> &r_mesh,
                                                               Real h,
                                                               int N,
                                                               Matrix<double> &density,
                                                               std::vector<Real> &tot_potential
) {

  tot_potential.resize(N);
  density.resize(N, 1);

  std::vector<Real> potential(N);
  std::vector<Real> potential_prev(N);
  std::vector<Real> vhartree_potential(N);
  std::vector<Real> dvhartree_potential(N);

  Matrix<Real> xc_potential(N, 1);
  Matrix<Real> e_xc(N, 1);
  Matrix<Real> density_prev(N, 1);

  std::vector<int> n;
  std::vector<int> l;
  std::vector<int> spin;
  std::vector<int> kappa;
  std::vector<Real> occupation;

  lsms::DefaultMixerVector mixVector(N);

  std::unique_ptr<lsms::AbstractMixer> mixer = std::make_unique<lsms::BroydenMixer>(mixVector.total_size,
                                                                                    alpha,
                                                                                    BROYDEN_MAX_STEP,
                                                                                    BROYDEN_ITER_RESET,
                                                                                    BROYDEN_W0);

  lsms::States::relativistic_atomic_states(
      Z, n, l, spin, kappa, occupation);

  std::vector<Real> e_ks_eigenenergies(n.size(), 0.0);

  if (iprint > 1) {
    fmt::print("n:    {:4d}\n", fmt::join(n, " "));
    fmt::print("l:    {:4d}\n", fmt::join(l, " "));
    fmt::print("spin: {:4d}\n", fmt::join(spin, " "));
    fmt::print("occ:  {:4.2f}\n", fmt::join(occupation, " "));
  }

  // Generate starting potential
  lsms::generate_starting_potential(tot_potential,
                                    r_mesh, Z,
                                    N
  );

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
    energy = lsms::generate_density(
        r_mesh,
        h,
        N,
        tot_potential,
        Z,
        n,
        l,
        spin,
        kappa,
        occupation,
        e_ks_eigenenergies,
        density,
        max_eig_iter,
        e_eig_tol
    );


    // XC potential
    xc->evaluate(r_mesh, h, density, N, e_xc, xc_potential);

    // Hartree potential
    lsms::radial_poisson(vhartree_potential,
                         dvhartree_potential,
                         r_mesh, h,
                         density, N);

    // Total energy
    tot_energy = lsms::total_energy(r_mesh,
                                    density,
                                    vhartree_potential,
                                    e_xc,
                                    tot_potential,
                                    energy,
                                    Z,
                                    N
    );

    // Mix density
    for (auto ir = 0; ir < N; ir++) {
      mixVector.vec_new[ir] = xc_potential[ir] * 0.5 + vhartree_potential[ir];
    }

    Real vrms = 0.0;
    Real qrms = 0.0;
    for (auto ir = 0; ir < N; ir++) {
      vrms += (mixVector.vec_new[ir] - mixVector.vec_old[ir]) * (mixVector.vec_new[ir] - mixVector.vec_old[ir]);
      qrms += (density[ir] - density_prev[ir]) * (density[ir] - density_prev[ir]);
    }

    mixVector.rms = std::sqrt(vrms);
    qrms = std::sqrt(qrms);

    mixer->mix(mixVector);

    for (auto ir = 0; ir < N; ir++) {
      mixVector.vec_old[ir] = mixVector.vec_new[ir];
      density_prev[ir] = density[ir];
    }

    if (std::abs(tot_energy - tot_energy_prev) < e_tol) {
      if (iprint >= 1) {
        fmt::print("{:>10}: {:>17d}\n", "Z", Z);
        fmt::print("{:>10}: {:>17d}\n", "Iter", iter);
        fmt::print("{:>10}: {:17.10f}\n", "One ele", energy);
        fmt::print("{:>10}: {:17.10f}\n", "Energy", tot_energy);
        fmt::print("{:>10}: {:17.4e}\n", "VRMS", vrms);
        fmt::print("{:>10}: {:17.4e}\n", "QRMS", qrms);
      }
      break;
    }

    tot_energy_prev = tot_energy;

  }

  if (iter >= max_iter) {
    fmt::print("WARNING: atomic dft didn't converge!\n");
  }

  return {e_ks_eigenenergies, tot_energy};

}
lsms::AtomicDFT::AtomicDFT(std::vector<int> functional,
                           int max_iter,
                           int max_eig_iter,
                           Real e_tol,
                           Real e_eig_tol,
                           int iprint,
                           Real alpha) :
    max_iter{max_iter}, max_eig_iter{max_eig_iter}, e_tol{e_tol}, e_eig_tol{e_eig_tol}, iprint{iprint}, alpha{alpha} {

  if (functional[0] == 0) {
    xc = std::make_unique<lsms::XCLDA>(1, functional);
  }

  if (functional[0] == 1) {
    xc = std::make_unique<lsms::XCLibxc>(1, functional);
  }

}
