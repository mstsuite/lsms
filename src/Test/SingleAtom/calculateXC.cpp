/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
// calculate the exchange correlation potential for a charge density

#include "Real.h"
#include <cmath>
#include <vector>

// calculate the Wigner-Seitz desnity parameter r_s
// r_s = (4 pi rho / 3)^{-1/3}
Real densityParameterRs(Real rho)
{
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
  
Real exchangePotentialLDA(Real rho)
{
  // v_x = -
}


