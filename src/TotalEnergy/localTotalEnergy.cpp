
#include "localTotalEnergy.hpp"

#include <cmath>

#include "Misc/integrator.hpp"
#include "Misc/poisson.hpp"
#include "PhysicalConstants.hpp"

/**
 * Calculates the total energy `energy` for each site
 *
 * @param lsms system parameters
 * @param atom atom object
 * @param energy total energy
 * @param pressure total pressure
 */
void localTotalEnergy(LSMSSystemParameters &lsms, AtomData &atom, Real &energy,
                      Real &pressure) {
  Real eigenvalueSum = 0.0;
  Real kineticEnergy = 0.0;
  Real coulombEnergy = 0.0;
  Real xcEnergy = 0.0;
  Real ezpt = 0.0;
  Real tpzpt = 0.0;
  Real lsf_energy = 0.0;
  Real ks = 0.0;
  Real rSphere;
  std::vector<Real> integral(atom.r_mesh.size());
  std::vector<Real> integrand(atom.r_mesh.size());

  energy = 0.0;
  pressure = 0.0;

  int end;
  switch (lsms.mtasa) {
    case 1: {
      end = atom.jws;
      rSphere = atom.rmt;
      break;
    }
    case 2: {
      end = atom.jws;
      rSphere = atom.rmt;
      break;
    }
    default:
      end = atom.jmt;
      rSphere = atom.rmt;
  }

  /**
   * Calculate the zeropoint energy
   */
  zeropt_(&ezpt, &tpzpt, &atom.omegaWS, &atom.ztotss);

  // calculate kinetic energy T:
  // T = \sum_{Core} e_i                            -- (1)
  //   + \int^{E_F} e n(e) de                       -- (2)
  //   - \int \rho(r) (v_{Coulomb} + v_{xc}) d^3r   -- (3)
  //   + \int m(r) (B_{xc} + B_{external}) d^3      -- (4)

  /**
   * Kinetic energy contributions
   */
  if (lsms.n_spin_pola == 1) {
    eigenvalueSum = atom.ecorv[0] + atom.esemv[0];

    kineticEnergy = atom.ecorv[0] + atom.esemv[0];  // (1)
    kineticEnergy += atom.evalsum[0];               // (2)

    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] = (atom.rhoNew(i, 0) * atom.vr(i, 0)) / (atom.r_mesh[i]);
    }

  } else {  // spin polarized
    eigenvalueSum =
        atom.ecorv[0] + atom.ecorv[1] + atom.esemv[0] + atom.esemv[1];

    kineticEnergy =
        atom.ecorv[0] + atom.ecorv[1] + atom.esemv[0] + atom.esemv[1];  // (1)
    kineticEnergy += atom.evalsum[0] + atom.evalsum[1];                 // (2)

    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] = (atom.rhoNew(i, 0) * atom.vr(i, 0) +
                      atom.rhoNew(i, 1) * atom.vr(i, 1)) /
                     (atom.r_mesh[i]);
    }
  }

  ks = lsms::radialIntegral(integrand, atom.r_mesh, end);
  kineticEnergy -= ks;

  if (lsms.global.iprint >= 0) {
    printf("evssum                      = %35.25lf Ry\n", eigenvalueSum);
    printf("ks                          = %35.25lf Ry\n", ks);
    printf("kinetic Energy              = %35.25lf Ry\n", kineticEnergy);
  }

  // calculate Coulomb Energy E_C:
  // E_C = \int \rho(r) v_{Coulomb} d^3r          -- (5)
  //     + E_{Madelung}                           -- (6) // external to this
  //     routine
  // break the Coulomb integral in two parts:
  //  \int \rho(r) v_{Coulomb} d^3r
  //      = \int \rho(r) \int^r rho(r')/r' dr' dr -- (5a)
  //      + \int rho(r) Z/r dr                    -- (5b)

  /**
   * Calculation of Coulomb contribution (5a)
   */
  std::vector<double> vhartreederiv(atom.r_mesh.size(), 0.0);
  std::vector<double> vhartree(atom.r_mesh.size(), 0.0);
  std::vector<double> density(atom.r_mesh.size(), 0.0);

  if (lsms.n_spin_pola == 1) {
    for (auto i = 0; i < atom.r_mesh.size(); i++) {
      density[i] = atom.rhoNew(i, 0);
    }
  } else {
    for (auto i = 0; i < atom.r_mesh.size(); i++) {
      density[i] = (atom.rhoNew(i, 0) + atom.rhoNew(i, 1));
    }
  }


  lsms::radial_poisson(vhartree, vhartreederiv, atom.r_mesh, atom.h, density,
                       end);

  if (lsms.n_spin_pola == 1) {
    for (auto i = 0; i < atom.r_mesh.size(); i++) {
      integral[i] = vhartree[i] * atom.rhoNew(i, 0);
    }
  } else {
    for (auto i = 0; i < atom.r_mesh.size(); i++) {
      integral[i] = vhartree[i] * (atom.rhoNew(i, 0) + atom.rhoNew(i, 1));
    }
  }

  double erho = lsms::radialIntegral(integral, atom.r_mesh, end);

  if (lsms.global.iprint >= 0) {
    printf("erho                        = %35.25lf Ry\n", erho);
  }

  /**
   * Calculation of the Coulomb contribution (5b)
   */
  if (lsms.n_spin_pola == 1) {
    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] = 2.0 * (atom.rhoNew(i, 0)) * atom.ztotss / (atom.r_mesh[i]);
    }
  } else {  // spin polarized
    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] = 2.0 * (atom.rhoNew(i, 0) + atom.rhoNew(i, 1)) *
                     atom.ztotss / (atom.r_mesh[i]);
    }
  }
  Real ezrho = -lsms::radialIntegral(integrand, atom.r_mesh, end);

  //  FILE *fp;
  //  fp = fopen("example.txt","w");
  //  for(int ir = 0; ir < atom.r_mesh.size(); ir++) {
  //    std::fprintf(fp, "%.30e %.30e\n", integrand[ir], atom.r_mesh[ir]);
  //  }
  //  std::fclose(fp);

  if (lsms.global.iprint >= 0) {
    printf("ezrho                       = %35.25lf Ry\n", ezrho);
  }

  coulombEnergy = erho + ezrho;  // (5)
  if (lsms.global.iprint >= 0)
    printf("Coulomb Energy              = %35.25lf Ry\n", coulombEnergy);

  /**
   * Exchange-Correlation energy                  -- (7)
   */
  if (lsms.n_spin_pola == 1) {
    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] = (atom.rhoNew(i, 0) * atom.exchangeCorrelationEnergy(i, 0));
    }
  } else {  // spin polarized
    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] = (atom.rhoNew(i, 0) * atom.exchangeCorrelationEnergy(i, 0) +
                      atom.rhoNew(i, 1) * atom.exchangeCorrelationEnergy(i, 1));
    }
  }
  xcEnergy = lsms::radialIntegral(integrand, atom.r_mesh, end);

  if (lsms.global.iprint >= 0) {
    printf("Exchange-Correlation Energy = %35.25lf Ry\n", xcEnergy);
    printf("ezpt                        = %35.25lf Ry\n", ezpt);
  }

  /**
   * Longitudinal spin fluctuations
   */
  if (lsms.n_spin_pola == 2) {
    auto mag_mom = atom.mvalws;
    lsf_energy += -convertKtoRydberg * lsms.temperature *
                  atom.lsf_functional.entropy(mag_mom);
  }

  if (lsms.global.iprint >= 0) {
    printf("LSF energy                  = %35.25lf Ry\n", lsf_energy);
  }

  energy += kineticEnergy + coulombEnergy + xcEnergy + ezpt + lsf_energy;
}

/**
 * Calculates the total energy `energy` for each site
 *
 * @param lsms system parameters
 * @param atom atom object
 * @param energy total energy class
 */
void localTotalEnergy(LSMSSystemParameters &lsms, AtomData &atom,
                      lsms::DFTEnergy &dft_energy) {
  Real xcEnergy = 0.0;
  Real ezpt = 0.0;
  Real tpzpt = 0.0;
  Real lsf_energy = 0.0;
  Real rSphere = 0.0;
  std::vector<Real> integral(atom.r_mesh.size());
  std::vector<Real> integrand(atom.r_mesh.size());

  int spin_channels;

  int end;
  switch (lsms.mtasa) {
    case 1: {
      end = atom.jws;
      rSphere = atom.rmt;
      break;
    }
    case 2: {
      end = atom.jws;
      rSphere = atom.rmt;
      break;
    }
    default:
      end = atom.jmt;
      rSphere = atom.rmt;
  }

  /**
   * Vibrational Zero-Point Energies (ZPE)
   */
  zeropt_(&ezpt, &tpzpt, &atom.omegaWS, &atom.ztotss);
  dft_energy.zero_point = ezpt;

  /**
   * Kinetic energy T:
   * T = \sum_{Core} e_i                            -- (1)
   *   + \int^{E_F} e n(e) de                       -- (2)
   *   - \int \rho(r) (v_{Coulomb} + v_{xc}) d^3r   -- (3)
   *   + \int m(r) (B_{xc} + B_{external}) d^3      -- (4)
   */
  if (lsms.n_spin_pola == 1) {
    spin_channels = 1;
  } else {  // spin polarized
    spin_channels = 2;
  }

  for (int spin = 0; spin < spin_channels; spin++) {
    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] +=
          (atom.rhoNew(i, spin) * atom.vr(i, spin)) / (atom.r_mesh[i]);
    }

    dft_energy.core_eigen += atom.ecorv[spin];
    dft_energy.semicore_eigen += atom.esemv[spin];
    dft_energy.one_ele += atom.evalsum[spin];
  }

  dft_energy.ks = lsms::radialIntegral(integrand, atom.r_mesh, end);

  dft_energy.kinetic = dft_energy.core_eigen + dft_energy.semicore_eigen +
                       dft_energy.one_ele - dft_energy.ks;

  /**
   * calculate Coulomb Energy E_C:
   * E_C = \int \rho(r) v_{Coulomb} d^3r          -- (5)
   *     + E_{Madelung}                           -- (6) // external to this
   *     routine
   * break the Coulomb integral in two parts:
   *  \int \rho(r) v_{Coulomb} d^3r
   *      = \int \rho(r) \int^r rho(r')/r' dr' dr -- (5a)
   *      + \int rho(r) Z/r dr                    -- (5b)
   *
   */

  std::vector<double> vhartreederiv(atom.r_mesh.size(), 0.0);
  std::vector<double> vhartree(atom.r_mesh.size(), 0.0);
  std::vector<double> density(atom.r_mesh.size(), 0.0);

  if (lsms.n_spin_pola == 1) {
    for (auto i = 0; i < atom.r_mesh.size(); i++) {
      density[i] = atom.rhoNew(i, 0);
    }
  } else {
    for (auto i = 0; i < atom.r_mesh.size(); i++) {
      density[i] = (atom.rhoNew(i, 0) + atom.rhoNew(i, 1));
    }
  }

  lsms::radial_poisson(vhartree, vhartreederiv, atom.r_mesh, atom.h, density,
                       end);

  if (lsms.n_spin_pola == 1) {
    for (auto i = 0; i < atom.r_mesh.size(); i++) {
      integral[i] = vhartree[i] * atom.rhoNew(i, 0);
    }
  } else {
    for (auto i = 0; i < atom.r_mesh.size(); i++) {
      integral[i] = vhartree[i] * (atom.rhoNew(i, 0) + atom.rhoNew(i, 1));
    }
  }

  double erho = lsms::radialIntegral(integral, atom.r_mesh, end);
  dft_energy.hartree = erho;

  /**
   * Calculation of the Coulomb contribution (5b)
   */
  if (lsms.n_spin_pola == 1) {
    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] = 2.0 * (atom.rhoNew(i, 0)) * atom.ztotss / (atom.r_mesh[i]);
    }
  } else {  // spin polarized
    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] = 2.0 * (atom.rhoNew(i, 0) + atom.rhoNew(i, 1)) *
                     atom.ztotss / (atom.r_mesh[i]);
    }
  }

  Real ezrho = -lsms::radialIntegral(integrand, atom.r_mesh, end);
  dft_energy.core_hartree = ezrho;

  dft_energy.coloumb = dft_energy.core_hartree + dft_energy.hartree;

  /**
   * Exchange-Correlation energy                  -- (7)
   */
  if (lsms.n_spin_pola == 1) {
    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] = (atom.rhoNew(i, 0) * atom.exchangeCorrelationEnergy(i, 0));
    }
  } else {  // spin polarized
    for (int i = 0; i < atom.r_mesh.size(); i++) {
      integrand[i] = (atom.rhoNew(i, 0) * atom.exchangeCorrelationEnergy(i, 0) +
                      atom.rhoNew(i, 1) * atom.exchangeCorrelationEnergy(i, 1));
    }
  }

  xcEnergy = lsms::radialIntegral(integrand, atom.r_mesh, end);
  dft_energy.xc = xcEnergy;

  /**
   * Longitudinal spin fluctuations
   */
  if (lsms.n_spin_pola == 2) {
    auto mag_mom = atom.mvalws;
    lsf_energy += -convertKtoRydberg * lsms.temperature *
                  atom.lsf_functional.entropy(mag_mom);
  }
  dft_energy.lsf = lsf_energy;

  dft_energy.total = dft_energy.kinetic + dft_energy.coloumb + dft_energy.xc +
                     dft_energy.zero_point + dft_energy.lsf;



}
