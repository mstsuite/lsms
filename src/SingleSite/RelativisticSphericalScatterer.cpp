/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include "Complex.hpp"
#include "SingleSiteScattering.hpp"
#include "Array3d.hpp"

#include "RelativisticScatterer.hpp"

#include "associatedLegendreFunction.hpp"

std::vector<int> RelativisticAngularMomentumIndices::kappaFromLambda;
std::vector<int> RelativisticAngularMomentumIndices::twoMuFromLambda;

void RelativisticAngularMomentumIndices::init(int lmax)
{
  int kmymax = 2*(lmax+1)*(lmax+1);
  kappaFromLambda.resize(kmymax);
  twoMuFromLambda.resize(kmymax);
  for(int kappa = -lmax - 1; kappa < lmax+1; kappa++)
  {
    if(kappa != 0)
    {
      int twoJ = 2*abs(kappa)-1;
      for(int twoMu = -twoJ; twoMu <= twoJ; twoMu += 2)
      {
        kappaFromLambda[lambdaIndexFromKappaMu(kappa, twoMu)] = kappa;
        twoMuFromLambda[lambdaIndexFromKappaMu(kappa, twoMu)] = twoMu;
      }
    }
  }
}

// lambdaIndexFromKappaMu
// kappa != 0
// twoMu = 2*mu, i.e. odd integers ..., -7, -5, -3, -1, +1, +3, +5, +7, ...
// inline int lambdaIndexFromKappaMu(int kappa, int twoMu)

// clebschGordonCoefficientHalf
// Clebsch-Gordon coefficients for j_2 = 1/2
// 
// ms: spin m quantum number. ms=+1/2 -> +1; ms=-1/2 -> -1
// inline double clebschGordonCoefficientJ2Half(int kappa, int m, int ms)


// clebschGordonCoefficientHalf
// Clebsch-Gordon coefficients for j_2 = 1/2
// C(l 1/2 j; mu - ms, ms)
// twoMu: 2*mu, i.e. ... -5, -3, -1, 1, 3, 5, ... odd integers
// ms: spin m quantum number. ms=+1/2 -> +1; ms=-1/2 -> -1
// inline double clebschGordonCoefficientJ2HalfTwoMu(int kappa, int twoMu, int ms)


// coupling of spin angular functions through sigma_z:
// int chi_kappa1,mu1^+ sigma_z chi_kappa2,mu2
// see P. Strange p. 427, eq. 11.80

// inline double spinAngularFunctionCouplingSigmaZ(int kappa1, int twoMu1, int kappa2, int twoMu2)

// class ChiSimgaZChiTable

std::vector<double> ChiSimgaZChiTable::GzKappa;
std::vector<int> ChiSimgaZChiTable::lambdaMinusKappaMinus1;
std::vector<double> ChiSimgaZChiTable::GzMinusKappaMinus1;

void ChiSimgaZChiTable::init(int lmax)
{
  GzKappa.resize(2*(lmax+1)*(lmax+1));
  lambdaMinusKappaMinus1.resize(2*(lmax+1)*(lmax+1));
  GzMinusKappaMinus1.resize(2*(lmax+1)*(lmax+1));
  
  for(int lambda=0; lambda<2*(lmax+1)*(lmax+1); lambda++)
  {
    int kappa = RelativisticAngularMomentumIndices::kappaFromLambda[lambda];
    int twoMu = RelativisticAngularMomentumIndices::twoMuFromLambda[lambda];
    GzKappa[lambda] = spinAngularFunctionCouplingSigmaZ(kappa, twoMu, kappa, twoMu);

    lambdaMinusKappaMinus1[lambda] = lambdaIndexFromKappaMu(-kappa-1, twoMu);
    // if(kappa == -1)
    if((kappa < 0) && (-2 * kappa - 1 == std::abs(twoMu))) 
    {
      lambdaMinusKappaMinus1[lambda] = -1;
      GzMinusKappaMinus1[lambda] = 0.0;
      // } else if
    } else {
      GzMinusKappaMinus1[lambda] = spinAngularFunctionCouplingSigmaZ(kappa, twoMu, -kappa-1, twoMu);
    }
  }
}

// spinAngularFunction \chi_\Lambda(\hat{r})

void spinAngularFunction(int lambda, Real theta, Real phi, Complex *chi) // , Real *plm)
{
  const Complex i(0, 1);
  int kappa = RelativisticAngularMomentumIndices::kappaFromLambda[lambda];
  int twoMu = RelativisticAngularMomentumIndices::twoMuFromLambda[lambda];
  int l = lFromKappa(kappa);
  std::vector<Real> plm((l+1)*(l+2)/2);
  associatedLegendreFunctionNormalized(std::cos(theta), l, &plm[0]); // Y_lm(\theta, \phi) = \bar{P}_{lm}(\cos \theta) e^{i m \phi}

  chi[0] = clebschGordonCoefficientJ2HalfTwoMu(kappa, twoMu, +1) * plm[plmIdx(l, (twoMu - 1)/2)] * std::exp(i * phi * Real((twoMu - 1)/2));
  chi[1] = clebschGordonCoefficientJ2HalfTwoMu(kappa, twoMu, -1) * plm[plmIdx(l, (twoMu + 1)/2)] * std::exp(i * phi * Real((twoMu + 1)/2));
}

/*
Input needed: momentum/energy, kappa, mu, r_mesh, v(r), b(r)
Output: t, C, S
 */

// N.B. kineticEnergy = relativisticEnergy - m_e c^2 (m_e = 1/2 in Rydberg units)

void freeElectronJ(Complex kineticEnergy, int kappa, int mu, std::vector<Real> &r_mesh, std::vector<DiracSpinor> &jKM)
{
  int l, lBar;
  if(kappa < 0)
  {
    l = -kappa - 1;
    lBar = -kappa;
  } else {
    l = kappa;
    lBar = kappa - 1;
  }
  int lmax = std::max(l, lBar);
  std::vector<Real> plm(((lmax+1)*(lmax+2))/2);

  // void associatedLegendreFunctionNormalized(R x, int lmax, R *Plm)
  // $ Y_lm(\theta, \phi) = \bar{P}_{lm}(\cos \theta) e^{i m \phi} $
  
}

void freeElectronN(Complex kineticEnergy, int kappa, int mu, std::vector<Real> &r_mesh, std::vector<DiracSpinor> &nKM)
{
  int l, lBar;
  if(kappa < 0)
  {
    l = -kappa - 1;
    lBar = -kappa;
  } else {
    l = kappa;
    lBar = kappa - 1;
  }
  int lmax = std::max(l, lBar);
  std::vector<Real> plm(((lmax+1)*(lmax+2))/2);
}
