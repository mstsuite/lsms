/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#include "newFunctionalInterface.hpp"

#include "Misc/rationalFit.hpp"
#include "Real.hpp"

extern "C" {
void get_rho_(Real *rho_in, Real *rho_out, int *stride, Real *r_mesh, int *jmt);
}

double alpha2_VWN(double rs, double dz, int ispin, double *xcEnergy) {
  const double a[3] = {-0.0337740, 0.06218140, 0.03109070};
  const double b[3] = {1.131070, 3.727440, 7.060420};
  const double c[3] = {13.00450, 12.93520, 18.05780};
  const double x0[3] = {-0.00475840, -0.104980, -0.32500};
  const double cst = 1.923661050;
  const double aip = 0.916330590;
  const double fnot = 1.709920950;
  const double bip = 0.259921050;
  const double for3 = 4.0 / 3.0;
  const double thrd = 1.0 / 3.0;

  double cx[3], q[3], bxx[3], tbq[3], tbxq[3], bb[3];

  for (int i = 0; i < 3; i++) {
    cx[i] = x0[i] * x0[i] + b[i] * x0[i] + c[i];
    q[i] = std::sqrt(4.0 * c[i] - b[i] * b[i]);
    bxx[i] = b[i] * x0[i] / cx[i];
    tbq[i] = 2.0 * b[i] / q[i];
    tbxq[i] = tbq[i] + 4.0 * x0[i] / q[i];
    bb[i] = 4.0 * b[i] * (1.0 - x0[i] * (b[i] + 2.0 * x0[i]) / cx[i]);
  }

  double zp1 = 1.0 + dz;
  double zm1 = 1.0 - dz;
  double xr = std::sqrt(rs);
  double pex = -aip / rs;
  double xrsq = rs;

  double g[3], dg[3];

  for (int i = 0; i < 3; i++) {
    double qi = q[i];
    double txb = 2.0 * xr + b[i];
    double fx = xrsq + xr * b[i] + c[i];
    double arct = std::atan2(qi, txb);
    double dxs = (xr - x0[i]) * (xr - x0[i]) / fx;
    g[i] = a[i] * (std::log(xrsq / fx) + tbq[i] * arct -
                   bxx[i] * (std::log(dxs) + tbxq[i] * arct));
    dg[i] =
        a[i] * (2.0 / xr - txb / fx - bxx[i] * (2.0 / (xr - x0[i]) - txb / fx) -
                bb[i] / (qi * qi + txb * txb));
  }

  double ecp = g[1];
  double zp3 = std::cbrt(zp1);
  double zm3 = std::cbrt(zm1);
  double zp3m3 = zp3 - zm3;
  //     part of last term in vx   eq(13)
  double fx1 = 0.5 * for3 * pex * zp3m3;
  double z4 = dz * dz * dz * dz;
  double fz = cst * (std::pow(zp1, for3) + std::pow(zm1, for3) - 2.0);
  double beta = fnot * (g[2] - g[1]) / g[0] - 1.0;
  double ec = ecp + fz * g[0] * (1.0 + z4 * beta) / fnot;
  double ex = pex * (1.0 + fz * bip);
  double f3ex = for3 * ex;

  //     echange-correlation energy
  *xcEnergy = ec + ex;

  //     exchange potential
  double vxx[2];
  vxx[0] = f3ex + fx1 * zm1;
  vxx[1] = f3ex - fx1 * zp1;
  //     correlation potential
  double vcc = ec - xr *
                        ((1.0 - z4 * fz) * dg[1] + z4 * fz * dg[2] +
                         (1.0 - z4) * fz * dg[0] / fnot) /
                        6.0;

  double facc =
      4.0 * g[0] *
      (dz * dz * dz * fz * beta + (1.0 + beta * z4) * zp3m3 / (6.0 * bip)) /
      fnot;

  //     exch-corr. potential for each spin as called in newpot
  if (ispin == 0) {
    return vcc + zm1 * facc + vxx[0];
  }

  return vcc - zp1 * facc + vxx[1];
}

double rsFromRho(double rho) { return std::cbrt(3.0 / (4.0 * M_PI * rho)); }

int NewFunctionalInterface::init(int nSpin, int *xcFunctional) {
  needGradients = needLaplacian = needKineticEnergyDensity = needExactExchange =
      false;
  numFunctionals = 0;
  if (xcFunctional[0] != 2) return 1;  // not a new lsms functional

  needGradients = false;
  needExactExchange = false;
  needLaplacian = false;
  needKineticEnergyDensity = false;

  for (int i = 1; i < numFunctionalIndices; i++) {
    if (xcFunctional[i] >= 0) {
      spinPolarized = false;
      if (nSpin > 1) spinPolarized = true;

      switch (xcFunctional[i]) {
        case 2:
          functional[numFunctionals] = XCFunctional::VoskoWilkNusair;
          break;

        default:
          printf(
              "Unknown Functional in new LSMS functional for functional %d!\n",
              xcFunctional[i]);
          exit(1);
      }

      numFunctionals++;
    }
  }
  return 0;
}

void NewFunctionalInterface::evaluate(std::vector<Real> &rMesh,
                                      Matrix<Real> &rhoIn, int jmt, int nSpin,
                                      Matrix<Real> &xcEnergyOut,
                                      Matrix<Real> &xcPotOut) {
  // note: rho in lsms is stored as 4*pi * r^2 * rho
  int nPts = jmt + 1;
  std::vector<Real> rho(nSpin * nPts);
  std::vector<Real> dRho(nSpin * nPts);
  std::vector<Real> sigma(
      (2 * nSpin - 1) *
      nPts);  // contracted gradient (1 entry/point for non polarized 3 for spin
              // polarized)(see libxc documentation)
  std::vector<Real> xcPot(nSpin * (jmt + 1)), xcEnergy(nPts);
  std::vector<Real> vSigma((2 * nSpin - 1) *
                           (jmt + 1));  // derivative with respect to contracted
                                        // gradient (see libxc documentation)

  if (nSpin == 1) {
    for (int ir = 0; ir < nPts; ir++) {
      Real r = rMesh[ir];  // (rMesh[ir] + 2.0*rMesh[ir+1])/3.0;
      rho[ir] = rhoIn(ir, 0) / (4.0 * M_PI * r * r);
      xcEnergyOut(ir, 0) = 0.0;
      xcPotOut(ir, 0) = 0.0;
      xcEnergyOut(ir, 1) = 0.0;
      xcPotOut(ir, 1) = 0.0;
    }
  } else {
    for (int ir = 0; ir < nPts; ir++) {
      rho[ir * 2] = rhoIn(ir, 0) / (4.0 * M_PI * rMesh[ir] * rMesh[ir]);
      xcEnergyOut(ir, 0) = 0.0;
      xcPotOut(ir, 0) = 0.0;
      rho[ir * 2 + 1] = rhoIn(ir, 1) / (4.0 * M_PI * rMesh[ir] * rMesh[ir]);
      xcEnergyOut(ir, 1) = 0.0;
      xcPotOut(ir, 1) = 0.0;
    }
  }
  if (needGradients) {
    // calculate the contracted gradients. Note that rho is spherically
    // symmetric: grad(rho) = e_r * (d rho / d r)
    if (nSpin > 1) {
      // spin polarized:
      // spin up
      calculateDerivative(&rMesh[0], &rho[0], &dRho[0], nPts, 2, 2);
      // spin down
      calculateDerivative(&rMesh[0], &rho[1], &dRho[1], nPts, 2, 2);
      for (int ir = 0; ir < nPts; ir++) {
        sigma[ir * 3] = dRho[ir * 2] * dRho[ir * 2];
        sigma[ir * 3 + 1] = dRho[ir * 2] * dRho[ir * 2 + 1];
        sigma[ir * 3 + 2] = dRho[ir * 2 + 1] * dRho[ir * 2 + 1];
      }
    } else {
      // non spin polarized
      calculateDerivative(&rMesh[0], &rho[0], &dRho[0], nPts);
      for (int ir = 0; ir < nPts; ir++) {
        sigma[ir] = dRho[ir] * dRho[ir];
      }
    }
  }
  for (int i = 0; i < numFunctionals; i++) {
    for (int j = 0; j < xcEnergy.size(); j++) xcEnergy[j] = 0.0;

    switch (functional[i]) {
      case XCFunctional::VoskoWilkNusair:
        for (int ir = 0; ir < jmt; ir++) {
          Real rhoTotal = rho[ir * nSpin];
          Real dz = 0.0;
          if (nSpin > 1) {
            rhoTotal += rho[ir * nSpin + 1];
            dz = (rho[ir * nSpin] - rho[ir * nSpin + 1]) / rhoTotal;
          }

          Real rs = rsFromRho(rhoTotal);
          xcPotOut(ir, 0) += alpha2_VWN(rs, dz, 0, &xcEnergy[ir]);
          if (nSpin > 1)
            xcPotOut(ir, 1) += alpha2_VWN(rs, dz, 1, &xcEnergy[ir]);
        }
        break;
      default:
        printf(
            "Unsuported Functional in new LSMS functional for functional!\n");
        exit(1);
    }
    if (nSpin == 1) {
      for (int ir = 0; ir < jmt; ir++) {
        xcEnergyOut(ir, 0) += xcEnergy[ir];
      }
    } else {
      for (int ir = 0; ir < jmt; ir++) {
        xcEnergyOut(ir, 0) += xcEnergy[ir];
        xcEnergyOut(ir, 1) += xcEnergy[ir];
      }
    }
  }
}

void NewFunctionalInterface::evaluateSingle(Real *rhoIn, int nSpin,
                                            Real *xcEnergyOut, Real *xcPotOut) {
  Real sigma[3];  // contracted gradient (see libxc documentation)
  Real xcPot[2], xcEnergy;
  Real vSigma[3];  // derivative with respect to contracted gradient (see libxc
                   // documentation)

  if (needGradients) {
    sigma[0] = 0.0;
    sigma[1] = 0.0;
    sigma[2] = 0.0;
  }

  Real rho = rhoIn[0];
  Real dz = 0.0;
  if (nSpin > 1) {
    rho += rhoIn[1];
    dz = (rhoIn[0] - rhoIn[1]) / rho;
  }

  Real rs = rsFromRho(rho);

  for (int i = 0; i < numFunctionals; i++) {
    xcEnergy = xcPot[0] = xcPot[1] = 0.0;
    switch (functional[i]) {
      case XCFunctional::VoskoWilkNusair:
        xcPot[0] = alpha2_VWN(rs, dz, 0, &xcEnergy);
        if (nSpin > 1) xcPot[1] = alpha2_VWN(rs, dz, 1, &xcEnergy);
        break;
      default:
        printf(
            "Unsuported Functional in new LSMS functional for functional!\n");
        exit(1);
    }

    *xcEnergyOut += xcEnergy;
    xcPotOut[0] += xcPot[0];
    if (nSpin > 1) {
      xcPotOut[1] += xcPot[1];
    }
  }
}
