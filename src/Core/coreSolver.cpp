/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#include "Real.hpp"
#include "Matrix.hpp"
#include "Main/SystemParameters.hpp"
#include "CoreStates.hpp"
#include "PhysicalConstants.hpp"

#include <vector>
#include <cmath>

extern "C"
{
  void deepst_(int *nqn, int *lqn, int *kqn, Real *en, Real *rv, Real *r, Real *rf, Real *h, Real *z, Real *c,
               int *nitmax, Real *tol, int *nws, int *nlast, int *iter, int *iprpts, int *ipdeq);
  void semcst_(int *nqn, int *lqn, int *kqn, Real *en,Real *rv, Real *r, Real *rf, Real *h, Real *z, Real *c,
               int *nitmax, Real *tol, int *nmt, int *nws, int *nlast, int *iter, int *iprpts,int *ipdeq);
  void newint_(int *nr, Real *r, Real *f, Real *g, int *ip0);
}

/*
    getcor_(&lsms.n_spin_pola,&lsms.mtasa,
            &local.atom[i].jmt,&local.atom[i].jws,&local.atom[i].r_mesh[0],&local.atom[i].h,&local.atom[i].xstart,
            &local.atom[i].vr(0,0),
            &local.atom[i].numc,&local.atom[i].nc(0,0),&local.atom[i].lc(0,0),&local.atom[i].kc(0,0),&local.atom[i].ec(0,0),
            &local.atom[i].ztotss,&local.atom[i].zsemss,&local.atom[i].zcorss,
            &local.atom[i].ecorv[0],&local.atom[i].esemv[0],&local.atom[i].corden(0,0),&local.atom[i].semcor(0,0),
            &lsms.nrelc,
            &local.atom[i].qcpsc_mt,&local.atom[i].qcpsc_ws,&local.atom[i].mcpsc_mt,&local.atom[i].mcpsc_ws,
            &local_iprpts,&local_ipcore,
            &lsms.global.iprint,lsms.global.istop,32);

      subroutine getcor(n_spin_pola,mtasa,
     >                  jmt,jws,r_mesh,h,xstart,vr,
     >                  numc,nc,lc,kc,ec,
     >                  ztotss,zsemss,zcorss,
     >                  ecorv,esemv,corden,semcor,
!ebot,etopcor,
     >                  nrelc,
     >                  qcpsc_mt,qcpsc_ws,mcpsc_mt,mcpsc_ws,
     >                  iprpts, ipcore,
     >                  iprint,istop)
*/

void getCoreStates(LSMSSystemParameters &lsms, AtomData &atom)
{
  /*
      ndeepz=(zcorss/n_spin_pola+.5d0)
      last=iprpts
      if(mtasa.eq.0) then
        last2=last
      else
        last2=jws
      endif
  */
  // Matrix<char> coreStateType(atom.numc,2);
  
  int numDeepStates = atom.zcorss;
  if(lsms.n_spin_pola == 2) numDeepStates = (atom.zcorss+1)/lsms.n_spin_pola;
  int last = atom.vr.l_dim();
  int last2 = last;
  int ndeep = 0;
  int nnorm;
  
  Real c = cphot * std::pow(10.0, lsms.nrelc);
  
  // if(lsms.mtasa != 0)
  //   last2 = atom.jws;
  int jmt = atom.jmt;
  if(lsms.mtasa==1) jmt = atom.jws;

  atom.corden = 0.0;
  atom.semcor = 0.0;
  atom.ecorv[0] = atom.ecorv[1] = 0.0;
  atom.esemv[0] = atom.esemv[1] = 0.0;
  atom.qcpsc_mt = atom.qcpsc_ws = 0.0;
  atom.mcpsc_mt = atom.mcpsc_ws = 0.0;
  atom.movedToValence[0] = atom.movedToValence[1] = 0;

  if(atom.numc <= 0) return;
  int local_iprpts=atom.vr.l_dim();
  std::vector<Real> f(local_iprpts+2);
  std::vector<Real> rtmp(local_iprpts+2);
  std::vector<Real> qmp(local_iprpts+2);
  
  for(int is=0; is<lsms.n_spin_pola; is++)
  {
    ndeep = 0;
    for(int ic=0; ic<atom.numc; ic++)
    {
      int nitmax = 50;
      int ipdeq = 5;
      int iter;
      Real tol = 1.0e-10;
      Real fac1 = (3-lsms.n_spin_pola) * std::abs(atom.kc(ic, is));
      ndeep += (3-lsms.n_spin_pola) * std::abs(atom.kc(ic, is));

      nnorm = last;

      deepst_(&atom.nc(ic, is), &atom.lc(ic, is), &atom.kc(ic, is), &atom.ec(ic, is),
              &atom.vr(0,is), &atom.r_mesh[0], &f[1], &atom.h, &atom.ztotss, &c,
              &nitmax, &tol, &atom.jws, &last, &iter, &local_iprpts, &ipdeq);
      
      if(ndeep <= numDeepStates)
      {
	nnorm = last;

/*
	deepst_(&atom.nc(ic, is), &atom.lc(ic, is), &atom.kc(ic, is), &atom.ec(ic, is),
                &atom.vr(0,is), &atom.r_mesh[0], &f[1], &atom.h, &atom.ztotss, &c,
                &nitmax, &tol, &atom.jws, &last, &iter, &local_iprpts, &ipdeq);
	/*
c        -------------------------------------------------------------
         call deepst(nc(i),lc(i),kc(i),ecore(i),
     >               rv,r,f(1),h,z,c,nitmax,tol,jws,last,iter,
     >               iprpts,ipdeq)
c        -------------------------------------------------------------
	 */
	atom.coreStateType(ic, is) = 'C'; // deep core state
        if(atom.ec(ic, is) >= lsms.energyContour.ebot || iter < 0)
          atom.coreStateType(ic, is) = 'V';
      } else {
	nnorm = last2;
        semcst_(&atom.nc(ic, is), &atom.lc(ic, is), &atom.kc(ic, is), &atom.ec(ic, is),
                &atom.vr(0,is), &atom.r_mesh[0], &f[1], &atom.h, &atom.ztotss, &c,
                &nitmax, &tol, &jmt, &atom.jws, &last2, &iter, &local_iprpts, &ipdeq);
	/*
c           ----------------------------------------------------------
            call semcst(nc(i),lc(i),kc(i),ecore(i),
     >                  rv,r,f(1),h,z,c,nitmax,tol,jmt,jws,last2,iter,
     >                  iprpts,ipdeq)
c           ----------------------------------------------------------
	 */
	atom.coreStateType(ic, is) = 'S'; // semi core state
	if(atom.ec(ic, is) >= lsms.energyContour.ebot || iter < 0)
	  atom.coreStateType(ic, is) = 'V'; // shallow core state -> move to valence
      }
        
      f[0] = 0.0;
      rtmp[0] = 0.0;
      for(int j = 1; j<=nnorm; j++)
      {
	rtmp[j] = std::sqrt(atom.r_mesh[j-1]);
	f[j] = f[j]/atom.r_mesh[j-1];
      }
      // call newint(nnorm+1,rtmp,f,qmp,3)
      int nnp1 = nnorm+1;
      int three = 3;
      newint_(&nnp1, &rtmp[0], &f[0], &qmp[0], &three);
      Real gnrm = 1.0 / (2.0 * qmp[nnorm-1]);
      for(int j=1; j<local_iprpts+1; j++)
        f[j] = f[j] * gnrm * atom.r_mesh[j-1];
      
      /*


! else
! if(nrelc.eq.0) then
!   kappa=kc(i)
! else
!   kappa=lc(i)
! endif
! call srcore(nrelc.eq.0,ecore(i),f(1),nc(i),kappa,potc,nr,
!    &        r,nitmax,tol,last2,1.d0)
! endif

c     ================================================================
c     normalize the wavefunctions
c     ================================================================
      f(0)=0.d0
      rtmp(0)=0.d0
! meis 4Sep20      do j=1,last2
      do j=1,nnorm
	rtmp(j)=sqrt(r(j))
	f(j)=f(j)/r(j)
      enddo
c     ----------------------------------------------------------------
! meis 4Sep20      call newint(last2+1,rtmp,f,qmp,3)
      call newint(nnorm+1,rtmp,f,qmp,3)
c     ----------------------------------------------------------------
      gnrm=1.d0/(two*qmp(last2))
      do j=1,last2
         f(j)=f(j)*gnrm*r(j)
      enddo    
      */
      if(atom.coreStateType(ic, is) == 'S')
      {
        for(int j=0; j<last2; j++)
          atom.semcor(j,is) += fac1 * f[j+1];
        atom.esemv[is] += fac1 * atom.ec(ic, is);
      } else if(atom.coreStateType(ic, is) == 'C') {
        for(int j=0; j<last; j++)
          atom.corden(j,is) += fac1 * f[j+1];
        atom.ecorv[is] += fac1 * atom.ec(ic, is);
      } else if(atom.coreStateType(ic, is) == 'V') {
        atom.movedToValence[is] += (3-lsms.n_spin_pola) * std::abs(atom.kc(ic, is));
      }
      /*
         if(ndeep.gt.ndeepz)then
            do j=1,last2
               semden(j)=semden(j) + fac1*f(j)
            enddo
            esemv=esemv+ecore(i)*fac1
	 else
            do j=1,last
! meis 4Sep20             do j=1,last2
               corden(j)= corden(j) + fac1*f(j)
            enddo
            ecorv=ecorv+ecore(i)*fac1
	 endif
      enddo
       */
    }
  }

  if(lsms.global.iprint >= 0)
  {
    for(int is=0; is<lsms.n_spin_pola; is++)
    {
      printf("getCoreStates: spin index #%d\n", is);
      printf("Eigenvalues: n    l    k         energy      type\n");
    
      for(int ic=0; ic<atom.numc; ic++)
      {
        printf("            %2d   %2d  %3d   %20.13f    %c\n",
               atom.nc(ic,is), atom.lc(ic,is), atom.kc(ic,is),
               atom.ec(ic,is), atom.coreStateType(ic,is));
      }
      printf("          Energy     core: %16.8f\n", atom.ecorv[is]);
      printf("          Energy semicore: %16.8f\n", atom.esemv[is]);
    }
    printf("\n");
  }

  Real qsemmt = 0.0;
  Real qsemws = 0.0;
  Real qcormt = 0.0;
  Real qcorws = 0.0;
  std::vector<Real> wrk1(local_iprpts+2);
  std::vector<Real> wrk2(local_iprpts+2);
  if(atom.numc > 0)
  {
    
    for(int is=0; is<lsms.n_spin_pola; is++)
    {
      rtmp[0] = 0.0;
      wrk2[0] = 0.0;
      for(int j=1; j<=last2; j++)
      {
        wrk2[j] = atom.semcor(j-1,is) / atom.r_mesh[j-1];
        rtmp[j] = std::sqrt(atom.r_mesh[j-1]);
      }
      int last2p1 = last2+1;
      int three = 3;
      newint_(&last2p1, &rtmp[0], &wrk2[0], &wrk1[0], &three);
      qsemmt += 2.0 * wrk1[jmt];
      qsemws += 2.0 * wrk1[last2-1];
      atom.mcpsc_mt += 2.0 * (1-2*is) * wrk1[jmt];
      atom.mcpsc_ws += 2.0 * (1-2*is) * wrk1[last2-1];

      wrk2[0] = 0.0;
      for(int j=1; j<=last2; j++)
      {
        wrk2[j] = atom.corden(j-1,is) / atom.r_mesh[j-1];
      }
      newint_(&last2p1, &rtmp[0], &wrk2[0], &wrk1[0], &three);
      qcormt += 2.0 * wrk1[jmt];
      qcorws += 2.0 * wrk1[last2-1];
      atom.mcpsc_mt += 2.0 * (1-2*is) * wrk1[jmt];
      atom.mcpsc_ws += 2.0 * (1-2*is) * wrk1[last2-1];
    }
  }

  Real qcorout = atom.zsemss + atom.zcorss - qsemmt - qcormt;
  atom.qcpsc_mt = qcormt + qsemmt;
  atom.qcpsc_ws = qcorws + qsemws;

  if(lsms.global.iprint >= 0)
  {
    printf("getCoreStates: Muffin-tin core charge        = %16.11f\n", qcormt);
    printf("               Muffin-tin semicore charge    = %16.11f\n", qsemmt);
    printf("               Muffin-tin core+semi moment   = %16.11f\n", atom.mcpsc_mt);
    printf("               Wigner-Seitz core+semi moment = %16.11f\n", atom.mcpsc_ws);
    printf("               Interstitial charge core      = %16.11f\n", qcorout);
  }
}

