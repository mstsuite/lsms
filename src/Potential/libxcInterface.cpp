/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include "Main/SystemParameters.hpp"
#include "Real.hpp"
#include "libxcInterface.hpp"
#include "Misc/rationalFit.hpp"

#ifdef USE_LIBXC
#include <xc.h>

const bool alwaysSpinPolarized=false;

extern "C"
{
  void get_rho_(Real *rho_in, Real *rho_out, int *stride,
               Real *r_mesh, int *jmt);
}

// LibxcInterface is a singleton class to provide an interface to libxc
// LibxcInterface::functional[numFunctionalIndices-1];
//  LibxcInterface::numFunctionals=0;
//  LibxcInterface::needGradients=false; // the functional needs gradients of the density (for GGAs)
//  LibxcInterface::needLaplacian=false; // need laplacians of the density (for MetaGGAs)
//  LibxcInterface::needKineticEnergyDensity=false; // for MetaGGAs
//  LibxcInterface::needExactExchange=false; // for Hybrid Functionals

int LibxcInterface::init(int nSpin, int *xcFunctional)
{
  needGradients=needLaplacian=needKineticEnergyDensity=needExactExchange=false;
  numFunctionals=0;
  if(xcFunctional[0]!=1) return 1; // not a libxc functional
  for(int i=1; i<numFunctionalIndices; i++)
  {
    if(xcFunctional[i]>=0)
    {
      int nspin=XC_UNPOLARIZED; if(nSpin>1 || alwaysSpinPolarized) nspin=XC_POLARIZED;
      if(xc_func_init(&functional[numFunctionals], xcFunctional[i], nspin)!=0) return 1;
      switch(functional[numFunctionals].info->family)
      {
      case XC_FAMILY_LDA: break;
      case XC_FAMILY_GGA: needGradients=true; break;
      case XC_FAMILY_HYB_GGA: needGradients=true; needExactExchange=true; break;
      case XC_FAMILY_MGGA: needGradients=true; needLaplacian=true; needKineticEnergyDensity=true; break;
      default: printf("Unknown Functional family in libxc for functional %d!\n",xcFunctional[i]); exit(1);
      }
      numFunctionals++;
    }
  }
  return 0;
}

void LibxcInterface::evaluate(std::vector<Real> &rMesh, Matrix<Real> &rhoIn, int jmt, int nSpin, Matrix<Real> &xcEnergyOut, Matrix<Real> &xcPotOut)
{
  int spinSize=nSpin;
  if (alwaysSpinPolarized) spinSize=2;
// note: rho in lsms is stored as 4*pi * r^2 * rho
  std::vector<Real> rho(spinSize*(jmt+1));
  std::vector<Real> dRho(spinSize*(jmt+1));
  std::vector<Real> sigma((2*spinSize-1)*(jmt+1)); // contracted gradient (1 entry/point for non polarized 3 for spin polarized)(see libxc documentation)
  std::vector<Real> xcPot(spinSize*(jmt+1)), xcEnergy(jmt+1);
  std::vector<Real> vSigma((2*spinSize-1)*(jmt+1)); // derivative with respect to contracted gradient (see libxc documentation)
  
  if(nSpin==1)
  {
    if (alwaysSpinPolarized)
    {
      for(int ir=0; ir<=jmt; ir++)
      {
        rho[ir*2  ] = 0.5*rhoIn(ir,0) / (4.0*M_PI*rMesh[ir]*rMesh[ir]);
        xcEnergyOut(ir,0)=0.0; xcPotOut(ir,0)=0.0;
        rho[ir*2+1] = 0.5*rhoIn(ir,0) / (4.0*M_PI*rMesh[ir]*rMesh[ir]);
        xcEnergyOut(ir,1)=0.0; xcPotOut(ir,1)=0.0;
      }
    } else {
      for(int ir=0; ir<=jmt; ir++)
      {
        Real r = rMesh[ir]; // (rMesh[ir] + 2.0*rMesh[ir+1])/3.0;
        rho[ir] = rhoIn(ir,0) / (4.0*M_PI*r*r);
        xcEnergyOut(ir,0)=0.0; xcPotOut(ir,0)=0.0;
        xcEnergyOut(ir,1)=0.0; xcPotOut(ir,1)=0.0;
      }
    }
    // int one=1;
    // get_rho_(&rhoIn(0,0), &rho[0], &one, &rMesh[0], &jmt);
  } else {
    for(int ir=0; ir<=jmt; ir++)
    {
      rho[ir*2  ] = rhoIn(ir,0) / (4.0*M_PI*rMesh[ir]*rMesh[ir]);
      xcEnergyOut(ir,0)=0.0; xcPotOut(ir,0)=0.0;
      rho[ir*2+1] = rhoIn(ir,1) / (4.0*M_PI*rMesh[ir]*rMesh[ir]);
      xcEnergyOut(ir,1)=0.0; xcPotOut(ir,1)=0.0;
    }
  }
  if(needGradients)
  {
// calculate the contracted gradients. Note that rho is spherically symmetric: grad(rho) = e_r * (d rho / d r)
    if(nSpin>1 || alwaysSpinPolarized)
    {
// spin polarized:
// spin up
      calculateDerivative(&rMesh[0], &rho[0], &dRho[0], jmt+1, 2, 2);
// spin down
      calculateDerivative(&rMesh[0], &rho[1], &dRho[1], jmt+1, 2, 2);
      for(int ir=0; ir<=jmt; ir++)
      {
	sigma[ir*3]=   dRho[ir*2]*dRho[ir*2];
	sigma[ir*3+1]= dRho[ir*2]*dRho[ir*2+1];
	sigma[ir*3+2]= dRho[ir*2+1]*dRho[ir*2+1];
      }
    } else {
      // non spin polarized
      calculateDerivative(&rMesh[0], &rho[0], &dRho[0], jmt+1);
      for(int ir=0; ir<=jmt; ir++)
      {
	sigma[ir] = dRho[ir]*dRho[ir];
      }
    }
  }
  for(int i=0; i<numFunctionals; i++)
  {
    for(int j=0; j<xcEnergy.size(); j++)
      xcEnergy[j] = 0.0;
    for(int j=0; j<xcPot.size(); j++)
      xcPot[j] = 0.0;
    // for(int ir=0; ir<jmt; ir++)
    // {
    switch(functional[i].info->family)
    {
    case XC_FAMILY_LDA: xc_lda_exc_vxc(&functional[i], jmt, &rho[0], &xcEnergy[0], &xcPot[0]); break;
    // case XC_FAMILY_LDA: xc_lda_exc_vxc(&functional[i], 1, &rho[ir*nSpin], &xcEnergy[ir], &xcPot[ir*nSpin]); break;
    case XC_FAMILY_GGA: xc_gga_exc_vxc(&functional[i], jmt, &rho[0], &sigma[0], &xcEnergy[0], &xcPot[0], &vSigma[0]); break;
    // case XC_FAMILY_GGA: xc_gga_exc_vxc(&functional[i], 1, &rho[ir*nSpin], &sigma[ir*3], &xcEnergy[ir], &xcPot[ir*(2*nSpin-1)],
    //                                    &vSigma[ir*(2*nSpin-1)]); break;
    default: printf("Unsuported Functional family in libxc for functional %d!\n",functional[i].info->number); exit(1);
    }
    if(nSpin == 1)
    {
      if (alwaysSpinPolarized)
      {
        for(int ir=0; ir<jmt; ir++)
        {
// libxc returns results in Hartree we need Rydberg as our energy units, so multiply by two
          xcEnergyOut(ir,0) += 2.0 * xcEnergy[ir];
          xcPotOut(ir,0)    += 2.0 * xcPot[ir*2    ];
          xcEnergyOut(ir,1) += 2.0 * xcEnergy[ir];
          xcPotOut(ir,1)    += 2.0 * xcPot[ir*2 + 1];
        }
      } else {
        for(int ir=0; ir<jmt; ir++)
        {
// libxc returns results in Hartree we need Rydberg as our energy units, so multiply by two
          xcEnergyOut(ir,0 ) += 2.0 * xcEnergy[ir];
          xcPotOut(ir,0)     += 2.0 * xcPot[ir];
        }
      }
    } else {
      for(int ir=0; ir<jmt; ir++)
      {
// libxc returns results in Hartree we need Rydberg as our energy units, so multiply by two
        xcEnergyOut(ir,0) += 2.0 * xcEnergy[ir];
        xcPotOut(ir,0)    += 2.0 * xcPot[ir*2    ];
        xcEnergyOut(ir,1) += 2.0 * xcEnergy[ir];
        xcPotOut(ir,1)    += 2.0 * xcPot[ir*2 + 1];
      }
    }
  }
}

void LibxcInterface::evaluateSingle(Real *rhoIn, int nSpin, Real *xcEnergyOut, Real *xcPotOut)
{
  Real sigma[3]; // contracted gradient (see libxc documentation)
  Real xcPot[2], xcEnergy;
  Real vSigma[3]; // derivative with respect to contracted gradient (see libxc documentation)

  if(needGradients)
  {
    sigma[0]= 0.0;
    sigma[1]= 0.0;
    sigma[2]= 0.0;
  }
  for(int i=0; i<numFunctionals; i++)
  {
    xcEnergy = xcPot[0] = xcPot[1] = 0.0;
    switch(functional[i].info->family)
    {
    case XC_FAMILY_LDA: xc_lda_exc_vxc(&functional[i], 1, &rhoIn[0], &xcEnergy, &xcPot[0]); break;
    case XC_FAMILY_GGA: xc_gga_exc_vxc(&functional[i], 1, &rhoIn[0], &sigma[0], &xcEnergy, &xcPot[0], &vSigma[0]); break;
    default: printf("Unsuported Functional family in libxc for functional %d!\n",functional[i].info->number); exit(1);
    }

// libxc returns results in Hartree? we need Rydberg as our energy units, so multiply by two
    *xcEnergyOut += 2.0*xcEnergy;  xcPotOut[0] += 2.0*xcPot[0];
    if(nSpin>1) { xcPotOut[1] += 2.0*xcPot[1]; }
  }
}

#else // we don't link with libxc


int LibxcInterface::init(int nSpin, int *xcFunctional)
{ printf("libxc is not linked with this version of LSMS!\n"); exit(1); }
void LibxcInterface::evaluate(std::vector<Real> &rMesh, Matrix<Real> &rhoIn, int jmt, int nSpin, std::vector<Real> &xcEnergyOut, Matrix<Real> &xcPotOut)
{ printf("libxc is not linked with this version of LSMS!\n"); exit(1); }
void LibxcInterface::evaluateSingle(Real *rhoIn, int nSpin, Real *xcEnergyOut, Real *xcPotOut)
{ printf("libxc is not linked with this version of LSMS!\n"); exit(1); }
#endif
