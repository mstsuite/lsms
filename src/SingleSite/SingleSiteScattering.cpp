/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include "SingleSiteScattering.hpp"



void calculateSingleScattererSolution(LSMSSystemParameters &lsms, AtomData &atom,
                                      Matrix<Real> &vr,
                                      Complex energy, Complex prel, Complex pnrel,
                                      NonRelativisticSingleScattererSolution &solution)
{
  int iprpts=atom.r_mesh.size();
  solution.energy=energy;
  // Real r_sph=atom.r_mesh[atom.jws];
  // if(lsms.mtasa==0) r_sph=atom.r_mesh[atom.jmt];
  Real r_sph = atom.rInscribed;
  if(lsms.mtasa > 0) r_sph = atom.rws;

  //YingWai's check
  //if(lsms.global.iprint>=0)
  //  printf("Inside calculateSingleScatterSolution. r_sph = %15.8f\n", r_sph);
  int writePhaseShift=0;
  
  if(lsms.lsmsMode == LSMSMode::phase_shifts) writePhaseShift = 1;

  if(lsms.n_spin_pola==1) // non spin polarized
  {
    single_scatterer_nonrel_(&lsms.nrelv, &lsms.clight, &atom.lmax, &atom.kkrsz,
                             &energy,&prel,&pnrel,
                             &vr(0,0),&atom.r_mesh[0],&atom.h,&atom.jmt,&atom.jws,
                             &solution.tmat_l(0,0,0),&solution.matom(0,0),
                             &solution.zlr(0,0,0),&solution.jlr(0,0,0),
                             &r_sph,
                             &iprpts,
                             &lsms.global.iprint,lsms.global.istop,&writePhaseShift,32);
    int kkrszsqr=atom.kkrsz*atom.kkrsz;
    int one=1;
    BLAS::zcopy_(&kkrszsqr,&solution.tmat_l(0,0,0),&one,&solution.tmat_g(0,0),&one);
  } else {
    for(int is=0; is<lsms.n_spin_pola; is++)
    {
      //YingWai's check
      //printf("Before single_scatterer_nonrel.\n");
      single_scatterer_nonrel_(&lsms.nrelv, &lsms.clight, &atom.lmax, &atom.kkrsz,
                               &energy,&prel,&pnrel,
                               &vr(0,is),&atom.r_mesh[0],&atom.h,&atom.jmt,&atom.jws,
                               &solution.tmat_l(0,0,is),&solution.matom(0,is),
                               &solution.zlr(0,0,is),&solution.jlr(0,0,is),
                               &r_sph,
                               &iprpts,
                               &lsms.global.iprint,lsms.global.istop,&writePhaseShift,32);
      //printf("After single_scatterer_nonrel.\n");
    }
    if(lsms.n_spin_cant>1)
    {
      trltog_(&atom.kkrsz,&atom.kkrsz,&atom.ubr[0],&atom.ubrd[0],
              &solution.tmat_l(0,0,0), &solution.tmat_l(0,0,1),&solution.tmat_g(0,0));
    } else {
      int kkrszsqr=atom.kkrsz*atom.kkrsz;
      int one=1;
      BLAS::zcopy_(&kkrszsqr,&solution.tmat_l(0,0,0),&one,&solution.tmat_g(0,0),&one);
      BLAS::zcopy_(&kkrszsqr,&solution.tmat_l(0,0,1),&one,&solution.tmat_g(0,atom.kkrsz),&one);
    }
  }
}

void calculateScatteringSolutions(LSMSSystemParameters &lsms, std::vector<AtomData> &atom,
                                  Complex energy, Complex prel, Complex pnrel,
                                  std::vector<NonRelativisticSingleScattererSolution> &solution)
{
  // if(atom.size()>solution.size()) solution.resize(atom.size());
  for(int i=0; i<atom.size(); i++)
    calculateSingleScattererSolution(lsms,atom[i],atom[i].vr,energy,prel,pnrel,solution[i]);
}

//=============== Relativistic ==============

void calculateSingleScattererSolution(LSMSSystemParameters &lsms, AtomData &atom,
                                      Matrix<Real> &vr,
                                      Complex energy,
                                      RelativisticSingleScattererSolution &solution)
{
   int iprpts=atom.r_mesh.size();
   solution.energy=energy;

// first we transform the potential from the up-down form to el-mag form
   std::vector<Real> vrr, brr;
   Matrix<Real> boprr;

   vrr.resize(iprpts);
   brr.resize(iprpts);
   boprr.resize(iprpts,2);


   // first we transform the potential from the up-down form to el-mag form
   for(int ir=0; ir<atom.r_mesh.size(); ir++)
   {
     vrr[ir] = 0.5 * (vr(ir,0) + vr(ir,1));
     brr[ir] = 0.5 * (vr(ir,0) - vr(ir,1));
     boprr(ir,0) = 0.0;
     boprr(ir,1) = 0.0;
   }

   Complex psq = energy + energy*energy/(lsms.clight*lsms.clight);

   int kmymax = 2*atom.kkrsz;
   int vacuumId = 1;
   int iflag=1;
   Real soscal =1.0;
   Real v0 = 0.0;
   int ir = atom.jws; // +1;
    if(lsms.mtasa==0) ir=atom.jmt+1;
   if(atom.ztotss<0.5) vacuumId=0;

   solution.gz.resize(iprpts, 2, kmymax);
   solution.gj.resize(iprpts, 2, kmymax);
   solution.fz.resize(iprpts, 2, kmymax);
   solution.fj.resize(iprpts, 2, kmymax);

   printf("iprpts = %d    ir = %d\n", iprpts, ir);
   // exit(1);
   
   single_scatterer_rel_(&energy, &psq, &atom.lmax, &kmymax,
                         &vacuumId, &v0,
                         &vrr[0], &brr[0], &boprr(0,0),
                         &atom.h, &ir, &atom.r_mesh[ir-1],
                         &solution.tmat_g(0,0),
                         &solution.gz(0,0,0), &solution.fz(0,0,0),
                         &solution.gj(0,0,0), &solution.fj(0,0,0),
                         &solution.nuz[0], &solution.indz(0,0), &iflag, &soscal,
                         &iprpts, &lsms.global.iprint,lsms.global.istop,32);

    for(int j1=0; j1<solution.gz.n_row(); j1++)
      for(int j2=0; j2<solution.gz.n_col(); j2++)
        for(int j3=0; j3<solution.gz.n_slice(); j3++)
        {
          printf("gz(%d,%d,%d) = %f %f\n",j1,j2,j3,solution.gz(j1,j2,j3).real(), solution.gz(j1,j2,j3).imag());
          printf("gj(%d,%d,%d) = %f %f\n",j1,j2,j3,solution.gj(j1,j2,j3).real(), solution.gj(j1,j2,j3).imag());
          printf("fz(%d,%d,%d) = %f %f\n",j1,j2,j3,solution.fz(j1,j2,j3).real(), solution.fz(j1,j2,j3).imag());
          printf("fj(%d,%d,%d) = %f %f\n",j1,j2,j3,solution.fj(j1,j2,j3).real(), solution.fj(j1,j2,j3).imag());
        }
   exit(0);
   // now we have tinv but we need t:
   Complex detl;
   gjinv_(&solution.tmat_g(0,0), &kmymax, &kmymax,&detl);
   // t now is in the local frame of refference, now we have to
   // transform it into the global frame:
   tripmt_(&atom.dmat(0,0),&solution.tmat_g(0,0),&atom.dmatp(0,0),&kmymax,&kmymax,&kmymax);
}
