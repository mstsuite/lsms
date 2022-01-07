/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include "Complex.hpp"
#include "SingleSiteScattering.hpp"
#include "Array3d.hpp"

#include "RelativisticScatterer.hpp"

class relativisticAngularMomentumIndices {
  std::vector<int> kappaFromLambda;
  std::vector<int> twoMuFromLambda;
public:
  void init(int lmax);
};

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

void spinAngularFunction()
{
}

/*
Input needed: momentum/energy, kappa, mu, r_mesh, v(r), b(r)
Output: t, C, S
 */

// N.B. kineticEnergy = relativisticEnergy - m_e c^2 (m_e = 1/2 in Rydberg units)

void freeElectronJ(Complex kineticEnergy, int kappa, int mu, std::vector<Real> &r_mesh)
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
}

void freeElectronN(Complex kineticEnergy, int kappa, int mu, std::vector<Real> &r_mesh)
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
}
