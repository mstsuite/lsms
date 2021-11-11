/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

// comparing lsms builtin xc functionals and libxc

#include <vector>
#include <cmath>
#include <cstdio>
#include <xc.h>
// #include "../../../opt/include/xc.h"

double alpha2_vbh(double rs, double dz, int ispin, double *xcEnergy)
{
  const double ccp = 0.0450;
  const double rp = 21.0;
  const double ccf = 0.02250;
  const double rf = 52.9166820;

  const double fm = std::pow(2.0, 4.0/3.0) - 2.0;

  double sp = 1.0 - double(ispin) * 2.0;

  double fdz = (std::pow(1.0 + dz, 4.0/3.0) + std::pow(1.0 - dz, 4.0/3.0) - 2.0)/fm;

  double ex = -0.916330/rs;

  double exf = ex * std::cbrt(2.0);
  double xp = rs/rp;
  double xf = rs/rf;

  double gp = (1.0 + xp*xp*xp) * std::log(1.0 + 1.0/xp)
              - xp*xp + xp/2.0 - 1.0/3.0;

  double gf = (1.0 + xf*xf*xf) * std::log(1.0 + 1.0/xf)
              - xf*xf + xf/2.0 - 1.0/3.0;

  double exc  = ex - ccp*gp;
  double excf = exf - ccf*gf;

  double dedz = (4.0/3.0) * (excf - exc)
                * (std::cbrt(1.0 + dz) - std::cbrt(1.0 - dz)) / fm;
  double gpp = 3.0 * xp*xp * std::log(1.0 + 1.0/xp) - 1.0/xp + 1.5 - 3.0 * xp;
  double gfp = 3.0 * xf*xf * std::log(1.0 + 1.0/xf) - 1.0/xf + 1.5 - 3.0 * xf;

  double depd = -ex/rs - ccp/rp * gpp;
  double defd = -exf/rs - ccf/rf * gfp;
  double decd = depd + (defd - depd) * fdz;

//     exchange-correlation energy
  *xcEnergy = exc + (excf - exc) * fdz;

//     exchange-correlation potential
  return exc + (excf - exc) * fdz - rs * decd/3.0 + sp * (1.0 - sp * dz) * dedz;
}

double alpha2_VWN(double rs, double dz, int ispin, double *xcEnergy)
{
  const double a[3] = {-0.0337740, 0.06218140, 0.03109070};
  const double b[3] = {1.131070, 3.727440, 7.060420};
  const double c[3] = {13.00450, 12.93520, 18.05780};
  const double x0[3] = {-0.00475840, -0.104980, -0.32500};
  const double cst = 1.923661050;
  const double aip = 0.916330590;
  const double fnot = 1.709920950;
  const double bip = 0.259921050;
  const double for3 = 4.0/3.0;
  const double thrd = 1.0/3.0;

  double cx[3], q[3], bxx[3], tbq[3], tbxq[3], bb[3];

  for(int i=0; i<3; i++)
  {
    cx[i] = x0[i]*x0[i] + b[i]*x0[i] + c[i];
    q[i] = std::sqrt(4.0d * c[i] - b[i]*b[i]);
    bxx[i] = b[i] * x0[i] / cx[i];
    tbq[i] = 2.0 * b[i] / q[i];
    tbxq[i] = tbq[i] + 4.0 * x0[i] / q[i];
    bb[i] = 4.0 * b[i] * (1.0 - x0[i] * (b[i] + 2.0 * x0[i]) / cx[i]);
  }

  double zp1 = 1.0 + dz;
  double zm1 = 1.0 - dz;
  double xr = std::sqrt(rs);
  double pex = -aip/rs;
  double xrsq = rs;

  double g[3], dg[3];

  for(int i=0; i<3; i++)
  {
    double qi = q[i];
    double txb = 2.0 * xr + b[i];
    double fx = xrsq + xr * b[i] + c[i];
    double arct= std::atan2(qi, txb);
    double dxs = (xr - x0[i]) * (xr - x0[i]) / fx;
    g[i] = a[i] * (std::log(xrsq / fx) + tbq[i] * arct
                    - bxx[i] * (std::log(dxs) + tbxq[i] * arct));
    dg[i] = a[i] * (2.0/xr - txb/fx 
                     - bxx[i] * (2.0/(xr - x0[i]) - txb/fx) 
                     - bb[i] / (qi*qi + txb*txb));
  }

  double ecp = g[1];
  double zp3 = std::cbrt(zp1);
  double zm3 = std::cbrt(zm1);
  double zp3m3 = zp3 - zm3;
//     part of last term in vx   eq(13)
  double fx1 = 0.5 * for3 * pex * zp3m3;
  double z4 = dz*dz*dz*dz;
  double fz= cst * (std::pow(zp1, for3)
                    + std::pow(zm1, for3) - 2.0);
  double beta = fnot * (g[2] - g[1]) / g[0] - 1.0;
  double ec = ecp + fz * g[0] * (1.0 + z4 * beta)/fnot;
  double ex = pex * (1.0 + fz * bip);
  double f3ex = for3 * ex;

//     echange-correlation energy
  *xcEnergy = ec + ex;

//     exchange potential
  double vxx[2];
  vxx[0] = f3ex + fx1 * zm1;
  vxx[1] = f3ex - fx1 * zp1;
//     correlation potential
  double vcc = ec - xr * ((1.0 - z4*fz) * dg[1] + z4 * fz * dg[2]
                + (1.0 - z4) * fz * dg[0] / fnot) / 6.0;

  double facc= 4.0 * g[0] * (dz*dz*dz * fz * beta 
                         + (1.0 + beta * z4) * zp3m3 / (6.0 * bip))/fnot;

//     exch-corr. potential for each spin as called in newpot
  if(ispin == 0)
  {
    return vcc + zm1 * facc + vxx[0];
  }

  return vcc - zp1 * facc + vxx[1];
}

double rsFromRho(double rho)
{
  return std::cbrt(3.0 / (4.0*M_PI * rho));
}

int main()
{
  int N = 100;
  double rhoMax = 200.0;

  std::vector<double> rho(N);
  std::vector<double> xcEnergy_legacy(N);
  std::vector<double> xcEnergy_libxc(N);
  std::vector<double> xEnergy_legacy(N);
  std::vector<double> xEnergy_libxc(N);

  const int libxc_lda_x = 1;
  const int libxc_lda_vBH = 17;
  const int libxc_lda_VWN = 7;

  xc_func_type xFunctional, cFunctional;
  if(xc_func_init(&xFunctional, libxc_lda_x, XC_UNPOLARIZED) != 0){
    fprintf(stderr, "Functional '%d' not found\n", libxc_lda_x);
    return 1;
  }
  if(xc_func_init(&cFunctional, libxc_lda_VWN, XC_UNPOLARIZED) != 0){
    fprintf(stderr, "Functional '%d' not found\n", libxc_lda_vBH);
    return 1;
  }


  for(int i = 0; i < N; i++)
  {
    rho[i] = ((double)(i+1)) * (rhoMax / ((double)(N)));
  }

  xc_lda_exc(&xFunctional, N, &rho[0], &xEnergy_libxc[0]);
  xc_lda_exc(&cFunctional, N, &rho[0], &xcEnergy_libxc[0]);

  for(int i = 0; i < N; i++)
  {
    double rs = rsFromRho(rho[i]);
    xEnergy_legacy[i] = -0.916330/rs;
    alpha2_VWN(rs, 0.0, 0, &xcEnergy_legacy[i]);

    xEnergy_libxc[i] *= 2.0;
    xcEnergy_libxc[i] *= 2.0;

    xcEnergy_libxc[i] += xEnergy_libxc[i];
  }

  for(int i = 0; i < N; i++)
  {
    printf("%g %g   %g %g   %g %g   %g\n",rho[i], rsFromRho(rho[i]),
           xEnergy_legacy[i], xcEnergy_legacy[i], xEnergy_libxc[i], xcEnergy_libxc[i],
           xcEnergy_legacy[i] - xcEnergy_libxc[i]);
  }

  return 0;
}
