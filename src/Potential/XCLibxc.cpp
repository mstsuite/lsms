//
// Created by F.Moitzi on 24.04.2022.
//

#include "XCLibxc.hpp"
#include "rationalFit.hpp"
#include "diff.hpp"

#include <stdexcept>

lsms::XCLibxc::XCLibxc(int nSpin, int xcFunctional[3])
    : XCBase(nSpin, xcFunctional) {

  setup(nSpin, xcFunctional, xcFunctional + 3);


}

lsms::XCLibxc::XCLibxc(int nSpin, std::vector<int> xcFunctional)
    : XCBase(nSpin, xcFunctional) {

  setup(nSpin, std::begin(xcFunctional), std::end(xcFunctional));


}

void lsms::XCLibxc::evaluate(std::vector<Real> &rMesh,
                             std::vector<Real> &drMesh,
                             const Matrix<Real> &rhoIn, int jmt,
                             Matrix<Real> &xcEnergyOut,
                             Matrix<Real> &xcPotOut) {
  /*
   * note: rho in lsms is stored as 4*pi * r^2 * rho
   */

  // Calculate spin channel size
  int nSigma;
  if (_nSpin == 1) {
    nSigma = 1;
  } else {
    nSigma = 3;
  }

  Matrix<double> g_xc(jmt, _nSpin);

  std::vector<Real> rho(_nSpin *jmt);

  std::vector<Real> rho_up(jmt);
  std::vector<Real> rho_down(jmt);

  std::vector<Real> drho_up(jmt);
  std::vector<Real> drho_down(jmt);

  std::vector<Real> xcEnergy(jmt);

  std::vector<Real> sigma(nSigma * jmt);
  std::vector<Real> xcPot(nSigma * jmt);
  std::vector<Real> vSigma(nSigma * jmt);

  xcEnergyOut = 0.0;
  xcPotOut = 0.0;
  g_xc = 0.0;

  if (_nSpin == 1) {
    // Non spin-polarized case

    for (int ir = 0; ir < jmt; ir++) {
      rho[ir] = rhoIn(ir, 0) / (4.0 * M_PI * rMesh[ir] * rMesh[ir]);
    }

  } else {
    // Spin-polarized case

    for (int ir = 0; ir < jmt; ir++) {
      rho[ir * 2] = rhoIn(ir, 0) / (4.0 * M_PI * rMesh[ir] * rMesh[ir]);
      rho_up[ir] = rho[ir * 2];
      rho[ir * 2 + 1] = rhoIn(ir, 1) / (4.0 * M_PI * rMesh[ir] * rMesh[ir]);
      rho_down[ir] = rho[ir * 2 + 1];
    }
  }

  if (needGradients) {
    if (_nSpin == 1) {

      drho_up = lsms::derivative<double>(rho.data(), jmt);

      for (int ir = 0; ir <= jmt; ir++) {
        sigma[ir] = drho_up[ir] * drho_up[ir] / (drMesh[ir] * drMesh[ir]);
      }

    } else {

      drho_up = lsms::derivative<double>(rho_up.data(), jmt);
      drho_down = lsms::derivative<double>(rho_down.data(), jmt);

      for (int ir = 0; ir < jmt; ir++) {
        sigma[ir * 3] = drho_up[ir] * drho_up[ir] / (drMesh[ir] * drMesh[ir]);
        sigma[ir * 3 + 1] = drho_up[ir] * drho_down[ir] / (drMesh[ir] * drMesh[ir]);
        sigma[ir * 3 + 2] = drho_down[ir] * drho_down[ir] / (drMesh[ir] * drMesh[ir]);
      }

    }
  }

  for (int i = 0; i < numFunctionals; i++) {

    switch (functionals[i].get_functional().info->family) {
      case XC_FAMILY_LDA:
        xc_lda_exc_vxc(&functionals[i].get_functional(), jmt, rho.data(),
                       xcEnergy.data(), xcPot.data());
        break;
      case XC_FAMILY_GGA:
        xc_gga_exc_vxc(&functionals[i].get_functional(), jmt, rho.data(),
                       sigma.data(), xcEnergy.data(), xcPot.data(),
                       vSigma.data());
        break;
    }


    if (needGradients) {

      if (_nSpin == 1) {
        // non-spin polarized

        auto isp = 0;
        for (int ir = 0; ir < jmt; ir++) {
          g_xc(ir, isp) += -2.0 * vSigma[ir] * drho_up[ir] / drMesh[ir];
        }

      } else {

        int isp;

        // Spin-up channel
        isp = 0;
        for (int ir = 0; ir < jmt; ir++) {

          g_xc(ir, isp) +=
              -2.0 * vSigma[ir * 3] * drho_up[ir] / drMesh[ir]
              - vSigma[ir * 3 + 1] * drho_down[ir] / drMesh[ir];

        }

        // Spin-down channel
        isp = 1;
        for (int ir = 0; ir < jmt; ir++) {

          g_xc(ir, isp) +=
              -2.0 * vSigma[ir * 3 + 2] * drho_down[ir] / drMesh[ir]
              - vSigma[ir * 3 + 1] * drho_up[ir] / drMesh[ir];

        }


      }


    }

    /*
     * LDA Part
     */

    //  Conversion from Hartree to Rydberg
    if (_nSpin == 1) {
      // Non spin-polarized case

      for (int ir = 0; ir < jmt; ir++) {
        xcEnergyOut(ir, 0) += 2.0 * xcEnergy[ir];
        xcPotOut(ir, 0) += 2.0 * xcPot[ir];
      }

    } else {
      for (int ir = 0; ir < jmt; ir++) {
        xcEnergyOut(ir, 0) += 2.0 * xcEnergy[ir];
        xcPotOut(ir, 0) += 2.0 * xcPot[ir * 2];

        xcEnergyOut(ir, 1) += 2.0 * xcEnergy[ir];
        xcPotOut(ir, 1) += 2.0 * xcPot[ir * 2 + 1];
      }
    }
  }

  if (needGradients) {
    /*
     * GGA Part
     */


    if (_nSpin == 1) {

      std::vector<double> v_grad_corr(jmt);

      int isp;

      isp = 0;
      v_grad_corr = lsms::derivative<double>(&g_xc(0, isp), jmt);

      for (int ir = 0; ir < jmt; ir++) {
        xcPotOut(ir, isp) += 2.0 * (
            2.0 * g_xc(ir, isp) / rMesh[ir]
            + v_grad_corr[ir] / drMesh[ir]);
      }

    } else {

      std::vector<double> v_grad_corr;

      int isp;

      isp = 0;
      v_grad_corr = lsms::derivative<double>(&g_xc(0, isp), jmt);

      for (int ir = 0; ir < jmt; ir++) {
        xcPotOut(ir, isp) += 2.0 * (
            2.0 * g_xc(ir, isp) / rMesh[ir]
            + v_grad_corr[ir] / drMesh[ir]);
      }

      isp = 1;
      v_grad_corr = lsms::derivative<double>(&g_xc(0, isp), jmt);

      for (int ir = 0; ir < jmt; ir++) {
        xcPotOut(ir, isp) += 2.0 * (
            2.0 * g_xc(ir, isp) / rMesh[ir]
            + v_grad_corr[ir]);
      }

    }


  }


}

void lsms::XCLibxc::evaluate(const Real rhoIn[2], Real &xcEnergyOut,
                             Real xcPotOut[2]) {
  /*
   * note: rho in lsms is stored as 4*pi * r^2 * rho
   */

  // Calculate spin channel size
  int nSigma;
  if (_nSpin == 1) {
    nSigma = 1;
  } else {
    nSigma = 3;
  }

  Real xcEnergy;
  std::vector<Real> sigma(nSigma, 0.0);
  std::vector<Real> vSigma(nSigma);
  std::vector<Real> xcPot(_nSpin);

  xcEnergyOut = 0.0;

  if (_nSpin == 1) {
    xcPotOut[0] = 0.0;
  } else {
    xcPotOut[0] = 0.0;
    xcPotOut[1] = 0.0;
  }

  if (needGradients) {
    // Calculate gradient
  }

  for (int i = 0; i < numFunctionals; i++) {
    switch (functionals[i].get_functional().info->family) {
      case XC_FAMILY_LDA:
        xc_lda_exc_vxc(&functionals[i].get_functional(), 1, rhoIn, &xcEnergy,
                       xcPot.data());
        break;
      case XC_FAMILY_GGA:
        xc_gga_exc_vxc(&functionals[i].get_functional(), 1, rhoIn, sigma.data(),
                       &xcEnergy, xcPot.data(), vSigma.data());
        break;
    }

    //  Conversion from Hartree to Rydberg
    if (_nSpin == 1) {
      xcEnergyOut += 2.0 * xcEnergy;
      xcPotOut[0] += 2.0 * xcPot[0];
    } else {
      xcEnergyOut += 2.0 * xcEnergy;
      xcPotOut[0] += 2.0 * xcPot[0];
      xcPotOut[1] += 2.0 * xcPot[1];
    }
  }
}

const std::vector<lsms::XCFuncType> &lsms::XCLibxc::get_functionals() const {
  return functionals;
}

std::string lsms::XCLibxc::get_name() {

  std::stringstream ss;

  for (auto &xc: functionals) {
    ss << xc.get_functional().info->name;

    if (&xc != &functionals.back()) { ss << "+"; };

  }

  ss << " (libxc)";

  return ss.str();
}

const xc_func_type &lsms::XCFuncType::get_functional() const {
  return _func_type;
}
