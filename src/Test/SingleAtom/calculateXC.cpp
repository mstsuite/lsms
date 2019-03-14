/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
// calculate the exchange correlation potential for a charge density

#include "Real.hpp"
#include <cmath>
#include <vector>

// calculate the Wigner-Seitz desnity parameter r_s
// r_s = (4 pi rho / 3)^{-1/3}
Real densityParameterRs(Real rho)
{
  return std::pow( (4.0*M_PI/3.0) * rho, -1.0/3.0);
}

/*
      function alpha2(rs,dz,sp,iexch,exchg)

c     von barth-hedin  exch-corr potential
c     j. phys. c5,1629(1972)
c
      data     ccp,rp,ccf,rf/0.0450d+00,21.0d+00,
     >                       0.02250d+00,52.9166820d+00/
c
c
  10  continue
      fm=2.0d+00**(4.0d+00/3.0d+00)-2.0d+00
      fdz = ((1.0d+00+dz)**(4.0d+00/3.0d+00)
     >     +(1.0d+00-dz)**(4.0d+00/3.0d+00)-2.0d+00)/fm
      ex=-0.916330d+00/rs
      exf=ex*2.0d+00**0.333333330d+00
      xp=rs/rp
      xf=rs/rf
      gp = (1.0d+00+xp**3)*log(1.0d+00+1.0d+00/xp)
     >    -xp*xp +xp/2.0d+00 - 0.333333330d+00
      gf = (1.0d+00+xf**3)*log(1.0d+00+1.0d+00/xf)
     >    -xf*xf +xf/2.0d+00 - 0.333333330d+00
      exc = ex-ccp*gp
      excf=exf-ccf*gf
      dedz= (4.0d+00/3.0d+00)*(excf-exc)
     >     *((1.0d+00+dz)**(1.0d+00/3.0d+00)
     >     -(1.0d+00-dz)**(1.0d+00/3.0d+00))/fm
      gpp = 3.0d+00*xp*xp*log(1.0d+00+1.0d+00/xp)-1.0d+00/xp
     >     +1.50d+00-3.0d+00*xp
      gfp = 3.0d+00*xf*xf*log(1.0d+00+1.0d+00/xf)-1.0d+00/xf
     >     +1.50d+00-3.0d+00*xf
      depd=-ex/rs-ccp/rp*gpp
      defd=-exf/rs-ccf/rf*gfp
      decd=depd+(defd-depd)*fdz
c     exchange-correlation energy
      exchg= exc + (excf-exc)*fdz
c     exchange-correlation potential
      alpha2 = exc+(excf-exc)*fdz-rs*decd/3.0d+00
     >        +sp*(1.0d+00-sp*dz)*dedz
      return
c
*/
Real alpha2(Real rs, Real dz, Real sp, Real &eXC)
{
//    von barth-hedin  exch-corr potential
//    j. phys. c5,1629(1972)
  const Real ccp = 0.0450;
  const Real rp = 21.0;
  const Real ccf = 0.02250;
  const Real rf = 52.9166820;
  
  Real fm=std::pow(2.0, (4.0/3.0))-2.0;
  Real fdz = (std::pow((1.0+dz),(4.0/3.0))
              + std::pow((1.0-dz),(4.0d+00/3.0d+00)) - 2.0)/fm;
  Real ex = -0.916330/rs;
  Real exf = ex * std::pow(2.0,1.0/3.0);
  Real xp = rs/rp;
  Real xf = rs/rf;
  Real gp = std::pow(1.0+xp,3) * std::log(1.0 + 1.0/xp)
    - xp*xp + xp/2.0 - 0.333333330;
  Real gf = std::pow(1.0+xf,3) * std::log(1.0 + 1.0/xf)
    - xf*xf + xf/2.0 - 0.333333330;
  Real exc = ex - ccp*gp;
  Real excf = exf - ccf*gf;
  Real dedz= (4.0/3.0) * (excf-exc)
    * (std::pow(1.0 + dz, 1.0/3.0)
       - std::pow(1.0 - dz, 1.0/3.0))/fm;
  Real gpp = 3.0 * xp*xp * std::log(1.0 + 1.0/xp) - 1.0/xp
    + 1.50 - 3.0*xp;
  Real gfp = 3.0 * xf*xf * std::log(1.0 + 1.0/xf) - 1.0/xf
    + 1.50 - 3.0*xf;
  Real depd=-ex/rs-ccp/rp*gpp;
  Real defd=-exf/rs-ccf/rf*gfp;
  Real decd=depd+(defd-depd)*fdz;
//     exchange-correlation energy
  eXC = exc + (excf-exc)*fdz;
//     exchange-correlation potential
  return exc + (excf-exc) * fdz - rs*decd/3.0 + sp*(1.0 - sp*dz)*dedz;
}
    
Real exchangePotentialLDA(Real rho, Real r)
{
  Real eXC;
  return alpha2(std::pow(3.0*r*r/rho , 1.0/3.0), 0.0, 1.0, eXC);
}

void exchangePotentialLDA(std::vector<Real> &rho, std::vector<Real> &r_mesh, std::vector<Real> &vXC)
{
  for(int i=0; i<rho.size(); i++)
  {
    vXC[i] = exchangePotentialLDA(rho[i], r_mesh[i]);
  }
}

