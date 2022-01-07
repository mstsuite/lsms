/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_RELATIVISTIC_SCATTERER_HPP
#define LSMS_RELATIVISTIC_SCATTERER_HPP

#include "Complex.hpp"
#include "SingleSiteScattering.hpp"
#include "Array3d.hpp"


class RelativisticAngularMomentumIndices {
public:
  static std::vector<int> kappaFromLambda;
  static std::vector<int> twoMuFromLambda;
public:
  void init(int lmax);
};


inline int lFromKappa(int kappa)
{
  if(kappa < 0)
    return -kappa - 1;
  return kappa;
}

// lambdaIndexFromKappaMu
// kappa != 0
// twoMu = 2*mu, i.e. odd integers ..., -7, -5, -3, -1, +1, +3, +5, +7, ...
// inline int lambdaIndexFromKappaMu(int kappa, int twoMu)
// {
//   return 2*kappa*kappa + kappa + (twoMu - 1)/2;
// }

// clebschGordonCoefficientHalf
// Clebsch-Gordon coefficients for j_2 = 1/2
// 
// ms: spin m quantum number. ms=+1/2 -> +1; ms=-1/2 -> -1
// inline double clebschGordonCoefficientJ2Half(int kappa, int m, int ms)
// {
//   int l;
//   double c, twolp1;
//   if(kappa < 0)
//   {
//     l = -kappa - 1;
//     twolp1 = 2.0*l + 1.0;
//     if(ms < 0)
//     {
//       c = std::sqrt((l - m + 1.0) / twolp1);
//     } else {
//       c = std::sqrt((l + m + 1.0) / twolp1);
//     }
//   } else {
//     l = kappa;
//     twolp1 = 2.0*l + 1.0;
//     if(ms < 0)
//     {
//       c =  std::sqrt((double)(l + m) / twolp1);
//     } else {
//       c = -std::sqrt((double)(l - m) / twolp1);
//     }
//   }
//  
//   return c;
// }

// clebschGordonCoefficientHalf
// Clebsch-Gordon coefficients for j_2 = 1/2
// C(l 1/2 j; mu - ms, ms)
// twoMu: 2*mu, i.e. ... -5, -3, -1, 1, 3, 5, ... odd integers
// ms: spin m quantum number. ms=+1/2 -> +1; ms=-1/2 -> -1
inline double clebschGordonCoefficientJ2HalfTwoMu(int kappa, int twoMu, int ms)
{
  int l;
  double c, twolp1;
  if(kappa < 0)
  {
    l = -kappa - 1;
    twolp1 = 2.0*l + 1.0;
    if(ms < 0)
    {
      c = std::sqrt((double)(l - (twoMu - 1)/2.0) / twolp1);
    } else {
      c = std::sqrt((double)(l + (twoMu + 1)/2.0) / twolp1);
    }
  } else {
    l = kappa;
    twolp1 = 2.0*l + 1.0;
    if(ms < 0)
    {
      c =  std::sqrt((double)(l + (twoMu + 1)/2.0) / twolp1);
    } else {
      c = -std::sqrt((double)(l - (twoMu - 1)/2.0) / twolp1);
    }
  }
  
  return c;
}

// coupling of spin angular functions through sigma_z:
// int chi_kappa1,mu1^+ sigma_z chi_kappa2,mu2
// see P. Strange p. 427, eq. 11.80

inline double spinAngularFunctionCouplingSigmaZ(int kappa1, int twoMu1, int kappa2, int twoMu2)
{
  int l1 = lFromKappa(kappa1);
  int l2 = lFromKappa(kappa2);

  if(l1 != l2)
    return 0.0;

  if(twoMu1 != twoMu2)
    return 0.0;
  
  return clebschGordonCoefficientJ2HalfTwoMu(kappa1, twoMu1, 1) * clebschGordonCoefficientJ2HalfTwoMu(kappa2, twoMu2, 1)
    - clebschGordonCoefficientJ2HalfTwoMu(kappa1, twoMu1, -1) * clebschGordonCoefficientJ2HalfTwoMu(kappa2, twoMu2, -1);
}

class ChiSimgaZChiTable
{
public:
  static std::vector<double> GzKappa;
  static std::vector<int> lambdaMinusKappaMinus1;
  static std::vector<double> GzMinusKappaMinus1;

  static void init(int lmax);
};

// spinAngularFunction \chi_\Lambda(\hat{r})

void spinAngularFunction();

/*
Input needed: momentum/energy, kappa, mu, r_mesh, v(r), b(r)
Output: t, C, S
 */

// N.B. kineticEnergy = relativisticEnergy - m_e c^2 (m_e = 1/2 in Rydberg units)

void freeElectronJ(Complex kineticEnergy, int kappa, int mu, std::vector<Real> &r_mesh);
void freeElectronN(Complex kineticEnergy, int kappa, int mu, std::vector<Real> &r_mesh);

#endif
