/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#include "coreSolver.hpp"

#include <cmath>
#include <vector>

#include "fmt/core.h"
#include "fmt/printf.h"

#include "Main/SystemParameters.hpp"
#include "Matrix.hpp"
#include "PhysicalConstants.hpp"
#include "Real.hpp"
#include "integrator.hpp"

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

//       subroutine getcor(n_spin_pola,mtasa,
//      >                  jmt,jws,r_mesh,h,xstart,vr,
//      >                  numc,nc,lc,kc,ec,
//      >                  ztotss,zsemss,zcorss,
//      >                  ecorv,esemv,corden,semcor,
// !ebot,etopcor,
//      >                  nrelc,
//      >                  qcpsc_mt,qcpsc_ws,mcpsc_mt,mcpsc_ws,
//      >                  iprpts, ipcore,
//      >                  iprint,istop)
*/

void getCoreStates(LSMSSystemParameters &lsms, AtomData &atom) {
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
  if (lsms.n_spin_pola == 2)
    numDeepStates = (atom.zcorss + 1) / lsms.n_spin_pola;
  int last = atom.vr.l_dim();
  int last2 = last;
  int ndeep = 0;
  int nnorm;

  Real c = cphot * std::pow(10.0, lsms.nrelc);

  // if(lsms.mtasa != 0)
  //   last2 = atom.jws;
  int jmt = atom.jmt;
  if (lsms.mtasa == 1) jmt = atom.jws;

  atom.corden = 0.0;
  atom.semcor = 0.0;
  atom.ecorv[0] = atom.ecorv[1] = 0.0;
  atom.esemv[0] = atom.esemv[1] = 0.0;
  atom.qcpsc_mt = atom.qcpsc_ws = 0.0;
  atom.mcpsc_mt = atom.mcpsc_ws = 0.0;
  atom.movedToValence[0] = atom.movedToValence[1] = 0;

  if (atom.numc <= 0) return;
  int local_iprpts = atom.vr.l_dim();
  std::vector<Real> f(local_iprpts + 2);
  std::vector<Real> rtmp(local_iprpts + 2);
  std::vector<Real> qmp(local_iprpts + 2);

  for (int is = 0; is < lsms.n_spin_pola; is++) {
    ndeep = 0;
    for (int ic = 0; ic < atom.numc; ic++) {
      int nitmax = 50;
      int ipdeq = 5;
      int iter;
      Real tol = 1.0e-10;
      Real fac1 = (3 - lsms.n_spin_pola) * std::abs(atom.kc(ic, is));
      ndeep += (3 - lsms.n_spin_pola) * std::abs(atom.kc(ic, is));

      nnorm = last;

      deepst_(&atom.nc(ic, is), &atom.lc(ic, is), &atom.kc(ic, is),
              &atom.ec(ic, is), &atom.vr(0, is), &atom.r_mesh[0], &f[1],
              &atom.h, &atom.ztotss, &c, &nitmax, &tol, &atom.jws, &last, &iter,
              &local_iprpts, &ipdeq);

      if (ndeep <= numDeepStates) {
        nnorm = last;
        atom.coreStateType(ic, is) = 'C';  // deep core state
        if (atom.ec(ic, is) >= lsms.energyContour.ebot || iter < 0)
          atom.coreStateType(ic, is) = 'V';
      } else {
        nnorm = last2;
        semcst_(&atom.nc(ic, is), &atom.lc(ic, is), &atom.kc(ic, is),
                &atom.ec(ic, is), &atom.vr(0, is), &atom.r_mesh[0], &f[1],
                &atom.h, &atom.ztotss, &c, &nitmax, &tol, &jmt, &atom.jws,
                &last2, &iter, &local_iprpts, &ipdeq);

        // std::printf("S: [%3d %3d] :%3d %3d %3d %20.10f\n", ic, is,
        // atom.nc(ic, is), atom.lc(ic, is), atom.kc(ic, is), atom.ec(ic, is));

        atom.coreStateType(ic, is) = 'S';  // semi core state
        if (atom.ec(ic, is) >= lsms.energyContour.ebot || iter < 0)
          atom.coreStateType(ic, is) =
              'V';  // shallow core state -> move to valence
      }

      f[0] = 0.0;
      rtmp[0] = 0.0;
      for (int j = 1; j <= nnorm; j++) {
        rtmp[j] = std::sqrt(atom.r_mesh[j - 1]);
        f[j] = f[j] / atom.r_mesh[j - 1];
      }

      // call newint(nnorm+1,rtmp,f,qmp,3)
      int nnp1 = nnorm + 1;
      int three = 3;
      newint_(&nnp1, &rtmp[0], &f[0], &qmp[0], &three);
      Real gnrm = 1.0 / (2.0 * qmp[nnorm - 1]);
      for (int j = 1; j < local_iprpts + 1; j++)
        f[j] = f[j] * gnrm * atom.r_mesh[j - 1];

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
      if (atom.coreStateType(ic, is) == 'S') {
        for (int j = 0; j < last2; j++) atom.semcor(j, is) += fac1 * f[j + 1];
        atom.esemv[is] += fac1 * atom.ec(ic, is);
      } else if (atom.coreStateType(ic, is) == 'C') {
        for (int j = 0; j < last; j++) atom.corden(j, is) += fac1 * f[j + 1];
        atom.ecorv[is] += fac1 * atom.ec(ic, is);
      } else if (atom.coreStateType(ic, is) == 'V') {
        atom.movedToValence[is] +=
            (3 - lsms.n_spin_pola) * std::abs(atom.kc(ic, is));
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

  if (lsms.global.iprint >= 0) {
    for (int is = 0; is < lsms.n_spin_pola; is++) {
      printf("getCoreStates: spin index #%d\n", is);
      printf("Eigenvalues: n    l    k         energy      type\n");

      for (int ic = 0; ic < atom.numc; ic++) {
        printf("            %2d   %2d  %3d   %20.13f    %c\n", atom.nc(ic, is),
               atom.lc(ic, is), atom.kc(ic, is), atom.ec(ic, is),
               atom.coreStateType(ic, is));
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
  std::vector<Real> wrk1(local_iprpts + 2);
  std::vector<Real> wrk2(local_iprpts + 2);
  if (atom.numc > 0) {
    for (int is = 0; is < lsms.n_spin_pola; is++) {
      rtmp[0] = 0.0;
      wrk2[0] = 0.0;
      for (int j = 1; j <= last2; j++) {
        wrk2[j] = atom.semcor(j - 1, is) / atom.r_mesh[j - 1];
        rtmp[j] = std::sqrt(atom.r_mesh[j - 1]);
      }
      int last2p1 = last2 + 1;
      int three = 3;
      newint_(&last2p1, &rtmp[0], &wrk2[0], &wrk1[0], &three);
      qsemmt += 2.0 * wrk1[jmt];
      qsemws += 2.0 * wrk1[last2 - 1];
      atom.mcpsc_mt += 2.0 * (1 - 2 * is) * wrk1[jmt];
      atom.mcpsc_ws += 2.0 * (1 - 2 * is) * wrk1[last2 - 1];

      wrk2[0] = 0.0;
      for (int j = 1; j <= last2; j++) {
        wrk2[j] = atom.corden(j - 1, is) / atom.r_mesh[j - 1];
      }
      newint_(&last2p1, &rtmp[0], &wrk2[0], &wrk1[0], &three);
      qcormt += 2.0 * wrk1[jmt];
      qcorws += 2.0 * wrk1[last2 - 1];
      atom.mcpsc_mt += 2.0 * (1 - 2 * is) * wrk1[jmt];
      atom.mcpsc_ws += 2.0 * (1 - 2 * is) * wrk1[last2 - 1];
    }
  }

  Real qcorout = atom.zsemss + atom.zcorss - qsemmt - qcormt;
  atom.qcpsc_mt = qcormt + qsemmt;
  atom.qcpsc_ws = qcorws + qsemws;

  if (lsms.global.iprint >= 0) {
    printf("getCoreStates: Muffin-tin core charge        = %16.11f\n", qcormt);
    printf("               Muffin-tin semicore charge    = %16.11f\n", qsemmt);
    printf("               Muffin-tin core+semi moment   = %16.11f\n",
           atom.mcpsc_mt);
    printf("               Wigner-Seitz core+semi moment = %16.11f\n",
           atom.mcpsc_ws);
    printf("               Interstitial charge core      = %16.11f\n", qcorout);
  }
}

void lsms::getNewCoreStates(LSMSSystemParameters &lsms, AtomData &atom) {
  // Default
  atom.corden = 0.0;
  atom.semcor = 0.0;
  atom.ecorv[0] = atom.ecorv[1] = 0.0;
  atom.esemv[0] = atom.esemv[1] = 0.0;
  atom.qcpsc_mt = atom.qcpsc_ws = 0.0;
  atom.mcpsc_mt = atom.mcpsc_ws = 0.0;
  atom.movedToValence[0] = atom.movedToValence[1] = 0;

  int nDeepStates;
  int iter;
  int mult;

  int nSphere;
  Real rSphere;

  int nitmax = NITMAX;
  int ipdeq = IPDEQ;
  Real tol = TOL;

  int nMesh = atom.vr.l_dim();

  Real c = cphot * std::pow(10.0, lsms.nrelc);

  std::vector<Real> dens(nMesh);

  int totDeepStates = 0;

  if (lsms.n_spin_pola == 1) {
    totDeepStates = (int) atom.zcorss;
  } else if (lsms.n_spin_pola == 2) {
    totDeepStates = std::ceil(atom.zcorss / 2.0);
  }

  // Bounding sphere radius
  nSphere = atom.jws;
  rSphere = atom.r_mesh[nSphere - 1];

  for (int is = 0; is < lsms.n_spin_pola; is++) {
    nDeepStates = 0;

    for (int ic = 0; ic < atom.numc; ic++) {
      mult = (3 - lsms.n_spin_pola) * std::abs(atom.kc(ic, is));
      nDeepStates += mult;

      /**
       * Calculate core states
       */

      if (nDeepStates <= totDeepStates) {
        // Deep core state
        deepst_(&atom.nc(ic, is), &atom.lc(ic, is), &atom.kc(ic, is),
                &atom.ec(ic, is), &atom.vr(0, is), atom.r_mesh.data(),
                dens.data(), &atom.h, &atom.ztotss, &c, &nitmax, &tol,
                &atom.jws, &nMesh, &iter, &nMesh, &ipdeq);

        atom.coreStateType(ic, is) = 'C';  // deep core state

      } else {
        // First deep core value
        deepst_(&atom.nc(ic, is), &atom.lc(ic, is), &atom.kc(ic, is),
                &atom.ec(ic, is), &atom.vr(0, is), atom.r_mesh.data(),
                dens.data(), &atom.h, &atom.ztotss, &c, &nitmax, &tol,
                &atom.jws, &nMesh, &iter, &nMesh, &ipdeq);

        // Semi-core state
        semcst_(&atom.nc(ic, is), &atom.lc(ic, is), &atom.kc(ic, is),
                &atom.ec(ic, is), &atom.vr(0, is), atom.r_mesh.data(),
                dens.data(), &atom.h, &atom.ztotss, &c, &nitmax, &tol,
                &atom.jmt, &atom.jws, &nMesh, &iter, &nMesh, &ipdeq);

        atom.coreStateType(ic, is) = 'S';  // semi core state
      }

      if (atom.ec(ic, is) >= lsms.energyContour.ebot || iter < 0) {
        atom.coreStateType(ic, is) = 'V';
      }

      /**
       * Renormalize densities
       */

      // Normalize to unity
      Real norm = lsms::radialIntegral(dens, atom.r_mesh, nMesh);

      for (int ir = 0; ir < nMesh; ir++) {
        dens[ir] = dens[ir] / norm;
      }

      // Reajust to sphere radius
      norm = lsms::radialIntegral(dens, atom.r_mesh, nSphere);

      Real dq = 1 - norm;
      Real factor = dq * 3 / (rSphere * rSphere * rSphere);

      for (int ir = 0; ir < nSphere; ir++) {
        dens[ir] = dens[ir] + factor * atom.r_mesh[ir] * atom.r_mesh[ir];
      }

      // Set density to zero outside
      for (int ir = nSphere; ir < nMesh; ir++) {
        dens[ir] = 0.0;
      }

      /**
       * Add density
       */

      if (atom.coreStateType(ic, is) == 'S') {
        for (int ir = 0; ir < nMesh; ir++) {
          atom.semcor(ir, is) += mult * dens[ir];
        }

        atom.esemv[is] += mult * atom.ec(ic, is);
      } else if (atom.coreStateType(ic, is) == 'C') {
        for (int ir = 0; ir < nMesh; ir++) {
          atom.corden(ir, is) += mult * dens[ir];
        }

        atom.ecorv[is] += mult * atom.ec(ic, is);
      } else if (atom.coreStateType(ic, is) == 'V') {
        atom.movedToValence[is] +=
            (3 - lsms.n_spin_pola) * std::abs(atom.kc(ic, is));
      }
    }
  }

  if (lsms.global.iprint >= 0) {
    for (int is = 0; is < lsms.n_spin_pola; is++) {
      printf("getCoreStates: spin index #%d\n", is);
      printf("Eigenvalues: n    l    k         energy      type\n");

      for (int ic = 0; ic < atom.numc; ic++) {
        printf("            %2d   %2d  %3d   %20.13f    %c\n", atom.nc(ic, is),
               atom.lc(ic, is), atom.kc(ic, is), atom.ec(ic, is),
               atom.coreStateType(ic, is));
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

  if (atom.numc > 0) {
    for (int is = 0; is < lsms.n_spin_pola; is++) {
      for (int ir = 1; ir < nMesh; ir++) {
        dens[ir] = atom.semcor(ir, is) + atom.corden(ir, is);
      }

      /**
       * Semi+deep states
       */

      Real mtNorm = lsms::radialIntegral(dens, atom.r_mesh, atom.jmt);
      Real wsNorm = lsms::radialIntegral(dens, atom.r_mesh, atom.jws);

      if (lsms.global.debug_core_states) {
        fmt::print("  Qcore  [{}]: {:12.8f}\n", is+1, mtNorm);
      }

      qsemmt += 2.0 * mtNorm;
      qsemws += 2.0 * wsNorm;
      atom.mcpsc_mt += 2.0 * (1 - 2 * is) * mtNorm;
      atom.mcpsc_ws += 2.0 * (1 - 2 * is) * wsNorm;
    }
  }

  Real qcorout = atom.zsemss + atom.zcorss - qsemmt - qcormt;
  atom.qcpsc_mt = qcormt + qsemmt;
  atom.qcpsc_ws = qcorws + qsemws;

  if (lsms.global.iprint >= 0) {
    printf("getCoreStates: Muffin-tin core charge        = %16.11f\n", qcormt);
    printf("               Muffin-tin semicore charge    = %16.11f\n", qsemmt);
    printf("               Muffin-tin core+semi moment   = %16.11f\n",
           atom.mcpsc_mt);
    printf("               Wigner-Seitz core+semi moment = %16.11f\n",
           atom.mcpsc_ws);
    printf("               Interstitial charge core      = %16.11f\n", qcorout);
  }
}
