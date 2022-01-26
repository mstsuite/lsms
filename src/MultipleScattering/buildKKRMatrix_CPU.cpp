/* -*- mode: C++; c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#include "Complex.hpp"
#include "Matrix.hpp"
#include <vector>
#include <cmath>

#include "buildKKRMatrix.hpp"

#include "SingleSite/SingleSiteScattering.hpp"
#include "MultipleScattering.hpp"
#include "Misc/Indices.hpp"
#include "Misc/Coeficients.hpp"
#include "Misc/associatedLegendreFunction.hpp"
#include "Main/LSMSMode.hpp"

// we might want to distinguish between systems where all lmax (and consequently kkrsz_ns) are the same
// and systems with potential different lmax on different atoms and l steps

// #define COMPARE_ORIGINAL 1


inline void calculateHankel(Complex prel, Real r, int lend, Complex *hfn)
{
  const Complex sqrtm1(0.0, 1.0);
  Complex z=prel*r;
  hfn[0]=-sqrtm1;
  hfn[1]=-1.0-sqrtm1/z;
  for(int l=1; l<lend; l++)
  {
    hfn[l+1]=(2.0*l + 1.0) * hfn[l]/z - hfn[l-1];
  }
  /*
c             l+1
c     hfn = -i   *h (k*R  )*sqrt(E)
c                  l    ij
*/
  z=std::exp(sqrtm1*z)/r;
  for(int l=0; l<=lend; l++)
  {
    hfn[l] = -hfn[l] * z * IFactors::ilp1[l];     
  }
}

inline void calculateSinCosPowers(Real *rij, int lend, Real *sinmp, Real *cosmp)
{
  const Real ptol = 1.0e-6;
  Real pmag = std::sqrt(rij[0]*rij[0]+rij[1]*rij[1]);
  cosmp[0] = 1.0;
  sinmp[0] = 0.0;
  if(pmag>ptol)
  {
    cosmp[1] = rij[0]/pmag;
    sinmp[1] = rij[1]/pmag;
  } else {
    cosmp[1] = 0.0;
    sinmp[1] = 0.0;
  }
  for(int m=2; m<=lend; m++)
  {
    cosmp[m] = cosmp[m-1]*cosmp[1] - sinmp[m-1]*sinmp[1];
    sinmp[m] = sinmp[m-1]*cosmp[1] + cosmp[m-1]*sinmp[1];
  }
}

void setBGijCPU(LSMSSystemParameters &lsms, AtomData &atom, int ir1, int ir2, int iOffset, int jOffset, Matrix<Complex> &bgij)
{
  if(lsms.n_spin_cant == 1) return;

  int kkri=(atom.LIZlmax[ir1]+1)*(atom.LIZlmax[ir1]+1);
  int kkrj=(atom.LIZlmax[ir2]+1)*(atom.LIZlmax[ir2]+1);
  int kkrsz = atom.kkrsz;

  if(lsms.relativity != full)
  {
    for(int i=0; i<kkri; i++)
      for(int j=0; j<kkrj; j++)
      {
        bgij(iOffset + kkri + i, jOffset        + j) = 0.0; // bgij(iOffset + i, jOffset + j);
        bgij(iOffset        + i, jOffset + kkrj + j) = 0.0; // bgij(iOffset + i, jOffset + j);
        bgij(iOffset + kkri + i, jOffset + kkrj + j) = bgij(iOffset + i, jOffset + j);
      }
  } else {
    /*
            call relmtrx(gij,bgij,kkr1,kkr2)
            fac=psq/ce
            do i=1,kkr1_ns
              do j=1,kkr2_ns
                bgij(i,j)=fac*bgij(i,j)
              end do
            end do
    */
    printf("Fully relativistic calculation not yet implemented in 'MultipleScattering/buildKKRMatrix.cpp : setBGijCPU'\n");
    exit(1);
  }
}

void buildBGijCPU(LSMSSystemParameters &lsms, AtomData &atom, int ir1, int ir2, Real *rij,
                  Complex energy, Complex prel, int iOffset, int jOffset, Matrix<Complex> &bgij)
{
  Complex hfn[2*lsms.maxlmax + 1];
  Real sinmp[2*lsms.maxlmax + 1];
  Real cosmp[2*lsms.maxlmax + 1];
  // Real plm[((lsms.maxlmax+1) * (lsms.maxlmax+2)) / 2];
  Real plm[lsms.angularMomentumIndices.ndlm];
  Complex dlm[lsms.angularMomentumIndices.ndlj];
  Real r = std::sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
  int lmax1 = atom.LIZlmax[ir1];
  int lmax2 = atom.LIZlmax[ir2];
  int kkri=(lmax1+1)*(lmax1+1);
  int kkrj=(lmax2+1)*(lmax2+1);
  int lend = lmax1 + lmax2;

  Real pi4=4.0*2.0*std::asin(1.0);

  calculateHankel(prel, r, lend, hfn);

  Real cosTheta = rij[2]/r;
  associatedLegendreFunctionNormalized<Real>(cosTheta, lend, plm);
  // for associatedLegendreFunctionNormalized all clm[i] == 1.0
  // for(int j=0;j<ndlm_local;j++)
  //   plm[j]=clm[j]*plm[j];
  
  //     calculate cos(phi) and sin(phi) .................................
  // needs to be serial
  calculateSinCosPowers(rij, lend, sinmp, cosmp);

  // can be parallel
  int j=0;
  for(int l=0; l<=lend; l++)
  {
    int ll = l*(l+1);
    j = ll;
    ll = ll/2;
    Real m1m = 1.0;
    dlm[j] = hfn[l]*plm[ll];
    for(int m=1; m<=l; m++)
    {
      m1m = -m1m;
      Complex fac = plm[ll+m] * std::complex<Real>(cosmp[m],sinmp[m]);
      dlm[j-m] = hfn[l]*m1m*fac;
      dlm[j+m] = hfn[l]*std::conj(fac);
    }
  }
//     ================================================================
//     calculate g(R_ij)...............................................
  // for(int i=0; i<kkri*kkrj; i++) gij[i]=0.0;
  for(int i=0; i<kkri; i++)
    for(int j=0; j<kkrj; j++)
      bgij(iOffset + i, jOffset + j) = 0.0;
  
//     loop over l1,m1............................................
  for(int lm1=0; lm1<kkrj; lm1++)
  {
    int l1=AngularMomentumIndices::lofk[lm1];
    int m1=AngularMomentumIndices::mofk[lm1];
    
//        loop over l2,m2..............................................
    for(int lm2=0; lm2<kkri; lm2++)
    {
      int l2=AngularMomentumIndices::lofk[lm2];
      int m2=AngularMomentumIndices::mofk[lm2];
      /*
        ==========================================================
                            l2-l1
           illp(lm2,lm1) = i

           perform sum over l3 with gaunt # ......................
        ==========================================================
      */
      int m3=m2-m1;
      int llow=std::max(std::abs(m3),std::abs(l1-l2));
      if(std::abs(prel)==0.0) llow=l1+l2;
      for(int l3=l1+l2; l3>=llow; l3-=2)
      {
        int j=l3*(l3+1)+m3;
        // gij[lm2+lm1*kkri] = gij[lm2+lm1*kkri]+cgnt(l3/2,lm1,lm2)*dlm[j];
        bgij(iOffset + lm2, jOffset + lm1) += GauntCoeficients::cgnt(l3/2,lm1,lm2)*dlm[j];
      }
      // gij[lm2+lm1*kkri]=pi4*illp(lm2,lm1)*gij[lm2+lm1*kkri];
      bgij(iOffset + lm2, jOffset + lm1) *= pi4 * IFactors::illp(lm2,lm1);
    }
  }

#ifdef COMPARE_ORIGINAL
  int kkr1 = kkri;
  int kkr2 = kkrj;
  bool exitCompare = false;
  Matrix<Complex> gijTest(kkr1,kkr2);
  Matrix<Complex> bgijTest(2*kkr1, 2*kkr2);
  int lmax=lsms.maxlmax;
  int kkrsz=(lmax+1)*(lmax+1);
  makegij_(&atom.LIZlmax[ir1],&kkr1,&atom.LIZlmax[ir2],&kkr2,
           &lsms.maxlmax,&kkrsz,&lsms.angularMomentumIndices.ndlj,&lsms.angularMomentumIndices.ndlm,
           &prel,&rij[0],&sinmp[0],&cosmp[0],
           &sphericalHarmonicsCoeficients.clm[0],&plm[0],
           &gauntCoeficients.cgnt(0,0,0),&gauntCoeficients.lmax,
           &lsms.angularMomentumIndices.lofk[0],&lsms.angularMomentumIndices.mofk[0],
           &iFactors.ilp1[0],&iFactors.illp(0,0),
           &hfn[0],&dlm[0],&gijTest(0,0),
           &pi4,&lsms.global.iprint,lsms.global.istop,32);
  int idx=0;
  for(int i=0; i<kkri; i++)
    for(int j=0; j<kkrj; j++)
    {
      if(bgij(iOffset + i, jOffset + j) != gijTest(i,j))
      // if(bgij[idx] != gijTest[idx])
      {
        printf("buildBGijCPU [idx=%d]: bgij(%d + %d, %d + %d) [%g + %gi] != gijTest(%d, %d) [%g + %gi]\n", idx,
               iOffset, i, jOffset, j, bgij(iOffset + i, jOffset + j).real(), bgij(iOffset + i, jOffset + j).imag(),
               i, j, gijTest(i,j).real(), gijTest(i,j).imag());
        exitCompare = true;
      }
      idx++;
    }
  if(exitCompare) exit(1);
#endif
        
  setBGijCPU(lsms, atom, ir1, ir2, iOffset, jOffset, bgij);
  
#ifdef COMPARE_ORIGINAL
  Complex psq=prel*prel;
  int kkr1_ns = 2*kkr1;
  int kkr2_ns = 2*kkr2;
  int nrel_rel=0;
  if(lsms.relativity==full) nrel_rel=1;
  setgij_(&gijTest(0,0),&bgijTest(0,0),&kkr1,&kkr1_ns,&kkr2,&kkr2_ns,
          &lsms.n_spin_cant,&nrel_rel,&psq,&energy);
  idx=0;
  for(int i=0; i<2*kkri; i++)
    for(int j=0; j<2*kkrj; j++)
    {
      // if(bgij(iOffset + i, jOffset + j) != bgijTest(i,j))
      if(bgij[idx] != bgijTest[idx])
      {
        printf("buildBGijCPU  [idx=%d]: bgij(%d + %d, %d + %d) [%g + %gi] != bgijTest(%d, %d) [%g + %gi]\n", idx,
               iOffset, i, jOffset, j, bgij(iOffset + i, jOffset + j).real(), bgij(iOffset + i, jOffset + j).imag(),
               i, j, bgijTest(i,j).real(), bgijTest(i,j).imag());
        exitCompare = true;
      }
      idx++;
    }
  if(exitCompare) exit(1);

  if((ir1==1 && ir2==0) || (ir1==10 && ir2==0))
  {
    printf("ir1=%d, ir2=%d: bgij(0,0) = %g + %gi; bgijTest(0,0) = %g + %gi\n",
           ir1, ir2, bgij(0,0).real(), bgij(0,0).imag(), bgijTest(0,0).real(), bgijTest(0,0).imag());
    printf("    rij = %g %g %g;  prel=%g + %gi\n", rij[0],  rij[1], rij[2], prel.real(), prel.imag());
    printf("    kkr1 = %d; kkr2 = %d; kkrsz = %d\n", kkr1, kkr2, kkrsz);
  }

#endif
}


void buildKKRMatrixLMaxIdenticalCPU(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int ispin, int iie, Complex energy, Complex prel,
                                    Matrix<Complex> &m)
{
  int nrmat_ns = lsms.n_spin_cant*atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block

  Complex cmone = Complex(-1.0,0.0);
  Complex czero=0.0;

  Matrix<Complex> bgij(nrmat_ns, nrmat_ns);
  Matrix<Complex> bgijSmall(kkrsz_ns, kkrsz_ns);
  
  m = 0.0; bgij = 0.0;
  for(int i=0; i<nrmat_ns; i++) m(i,i)=1.0;

  // loop over the LIZ blocks
  for(int ir1 = 0; ir1 < atom.numLIZ; ir1++)
  {
    int iOffset = ir1 * kkrsz_ns; // this assumes that there are NO lStep reductions of lmax!!!
    for(int ir2 = 0; ir2 < atom.numLIZ; ir2++)
    {
      if(ir1 != ir2)
      {
        int jOffset = ir2 * kkrsz_ns; // this assumes that there are NO lStep reductions of lmax!!!
        Real rij[3];
        int lmax1 = atom.LIZlmax[ir1];
        int lmax2 = atom.LIZlmax[ir2];
        int kkr1=(lmax1+1)*(lmax1+1);
        int kkr2=(lmax2+1)*(lmax2+1);
        int kkr1_ns = kkr1 * lsms.n_spin_cant;
        int kkr2_ns = kkr2 * lsms.n_spin_cant;
        rij[0]=atom.LIZPos(0,ir1)-atom.LIZPos(0,ir2);
        rij[1]=atom.LIZPos(1,ir1)-atom.LIZPos(1,ir2);
        rij[2]=atom.LIZPos(2,ir1)-atom.LIZPos(2,ir2);
        
        buildBGijCPU(lsms, atom, ir1, ir2, rij, energy, prel, iOffset, jOffset, bgij);
        // buildBGijCPU(lsms, atom, ir1, ir2, rij, energy, prel, 0, 0, bgijSmall);

#ifdef COMPARE_ORIGINAL
        int kkri=(lmax1+1)*(lmax1+1);
        int kkrj=(lmax2+1)*(lmax2+1);
        int lmax=lsms.maxlmax;
        int kkrsz=(lmax+1)*(lmax+1);
        bool exitCompare = false;
        Matrix<Complex> tmat_n(lsms.n_spin_cant*atom.kkrsz, lsms.n_spin_cant*atom.kkrsz);
        int im=0;
        if(lsms.n_spin_pola == lsms.n_spin_cant) // non polarized or spin canted
        {
          for(int js=0; js<lsms.n_spin_cant; js++)
          {
            int jsm = kkrsz*kkrsz_ns*js;
            for(int j=0; j<kkr1; j++)
            {
              for(int is=0; is<lsms.n_spin_cant; is++)
              {
                int jm=jsm+kkrsz_ns*j+kkrsz*is;
                int one=1;
                BLAS::zcopy_(&kkr1,&local.tmatStore(iie*local.blkSizeTmatStore+jm,atom.LIZStoreIdx[ir1]),&one,&tmat_n[im],&one);
                im+=kkr1;
              }
            }
          }
        } else { // spin polarized colinear version for ispin
          int ispin=0;
          printf("warning: cant't test building kkrMatrix for collinear spin polarized yet!\n");
          exit(1);
          int jsm = kkrsz*kkrsz*ispin; // copy spin up or down?
          for(int j=0; j<kkr1; j++)
          {
            int jm=jsm+kkrsz_ns*j;
            int one=1;
            BLAS::zcopy_(&kkr1,&local.tmatStore(iie*local.blkSizeTmatStore+jm,atom.LIZStoreIdx[ir1]),&one,&tmat_n[im],&one);
            im+=kkr1;
          }
        }

        int idx = 0;
        for(int i=0; i<2*kkri; i++)
          for(int j=0; j<2*kkrj; j++)
          {
            int jm = i + j*2*kkri;
            if(tmat_n(i, j) != local.tmatStore(iie*local.blkSizeTmatStore+jm,atom.LIZStoreIdx[ir1]))
            // if(tmat_n[idx] != local.tmatStore(iie*local.blkSizeTmatStore+idx,atom.LIZStoreIdx[ir1]))
            {
              printf("buildKKRMatrix...CPU [idx=%d]: tmat_n(%d, %d) [%g + %gi] != tmatStore(%d + %d, %d) [%g + %gi]\n", idx,
                     i, j, tmat_n(i, j).real(), tmat_n(i, j).imag(),
                     iie*local.blkSizeTmatStore, jm, atom.LIZStoreIdx[ir1],
                     local.tmatStore(iie*local.blkSizeTmatStore+jm,atom.LIZStoreIdx[ir1]).real(),
                     local.tmatStore(iie*local.blkSizeTmatStore+jm,atom.LIZStoreIdx[ir1]).imag());
              exitCompare = true;
            }
            idx++;
          }

        if(exitCompare) exit(1);
        
        Matrix<Complex> bgijTest(kkr1_ns, kkr2_ns);
        for(int i=0; i<kkr1_ns; i++)
          for(int j=0; j<kkr2_ns; j++)
            bgijTest(i,j) = bgij(iOffset + i, jOffset + j);
        /*
        BLAS::zgemm_("n", "n", &kkr1_ns, &kkr2_ns, &kkr1_ns, &cmone,
                     &tmat_n(0, 0), &kkr1_ns,
                     &bgijTest(0, 0), &kkr1_ns, &czero,
                     &m(iOffset, jOffset), &nrmat_ns);
        */
        if((ir1==1 && ir2==0) || (ir1==10 && ir2==0))
        {
          Complex p = -tmat_n(0,0)*bgijSmall(0,0);
          printf("ir1=%d, ir2=%d: bgijSmall(0,0) = %g +%g i; tmat_n(0,0) = %g + %gi; -product =  %g + %gi\n",
                 ir1, ir2, bgijSmall(0,0).real(), bgijSmall(0,0).imag(), tmat_n(0,0).real(), tmat_n(0,0).imag(), p.real(), p.imag());
        }
#endif
        if(lsms.n_spin_pola == lsms.n_spin_cant) // non polarized or spin canted
        {
          BLAS::zgemm_("n", "n", &kkr1_ns, &kkr2_ns, &kkr1_ns, &cmone,
                       &local.tmatStore(iie*local.blkSizeTmatStore, atom.LIZStoreIdx[ir1]), &kkr1_ns,
                       // &tmat_n(0, 0), &kkr1_ns,
                       &bgij(iOffset, jOffset), &nrmat_ns, &czero,
                       // &bgijSmall(0, 0), &kkrsz_ns, &czero,
                       &m(iOffset, jOffset), &nrmat_ns);
        } else {  // spin polarized, collinear
          int lmax=lsms.maxlmax;
          int kkrsz=(lmax+1)*(lmax+1);
          int spinOffset = kkrsz * kkrsz * ispin;
          BLAS::zgemm_("n", "n", &kkr1_ns, &kkr2_ns, &kkr1_ns, &cmone,
                       &local.tmatStore(iie*local.blkSizeTmatStore + spinOffset, atom.LIZStoreIdx[ir1]), &kkr1_ns,
                       // &tmat_n(0, 0), &kkr1_ns,
                       &bgij(iOffset, jOffset), &nrmat_ns, &czero,
                       // &bgijSmall(0, 0), &kkrsz_ns, &czero,
                       &m(iOffset, jOffset), &nrmat_ns);
        }
        /*
        for(int i=0; i<kkr1_ns; i++)
          for(int j=0; j<kkr2_ns; j++)
          {
            m(iOffset + i, jOffset + j) = 0.0;
            for(int k=0; k<kkr1_ns ; k++)
              m(iOffset + i, jOffset + j) -= tmat_n(i, k) * // local.tmatStore(iie*local.blkSizeTmatStore + , atom.LIZStoreIdx[ir1]) *
                // bgij(iOffset + k, jOffset + j);
                bgijSmall(k, j);
          }
        */
        
      }
    }
  }

#ifdef COMPARE_ORIGINAL
  bool exitCompare = false;
  // int ispin = 0;
  Matrix<Complex> mTest(nrmat_ns, nrmat_ns);
  buildKKRMatrix(lsms, local, atom, ispin, energy, prel, iie, mTest);
  int idx=0;
  for(int j=0; j<nrmat_ns; j++)
    for(int i=0; i<nrmat_ns; i++)
    {
      if(m(i, j) != mTest(i,j))
      // if(m[idx] != mTest[idx])
      {
        printf("buildKKRMatrix...CPU [idx=%d]: m(%d, %d) [%g + %gi] != mTest(%d, %d) [%g + %gi]\n", idx,
               i, j, m(i, j).real(), m(i, j).imag(),
               i, j, mTest(i, j).real(), mTest(i, j).imag());
        exitCompare = true;
        m(i, j) = mTest(i,j);
      }
      idx++;
    }
  
  if(exitCompare) exit(1);
#endif
}

void buildKKRMatrixLMaxDifferentCPU(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int ispin, int iie, Complex energy, Complex prel,
                                    Matrix<Complex> &m)
{
  int nrmat_ns = lsms.n_spin_cant*atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block

  const Complex cmone=-1.0;
  const Complex czero=0.0;

  Matrix<Complex> bgij(nrmat_ns, nrmat_ns);
  
  m = 0.0; bgij = 0.0;
  for(int i=0; i<nrmat_ns; i++) m(i,i)=1.0;

  std::vector<int> offsets(atom.numLIZ);
  offsets[0] = 0;
  for(int ir = 1; ir < atom.numLIZ; ir++)
    offsets[ir] = offsets[ir-1] + lsms.n_spin_cant * (atom.LIZlmax[ir-1]+1)*(atom.LIZlmax[ir-1]+1);
  
  // loop over the LIZ blocks
  for(int ir1 = 0; ir1 < atom.numLIZ; ir1++)
  {
    int iOffset = offsets[ir1];
    for(int ir2 = 0; ir2 < atom.numLIZ; ir2++)
    {
      if(ir1 != ir2)
      {
        int jOffset = offsets[ir2];
        int lmax1 = atom.LIZlmax[ir1];
        int lmax2 = atom.LIZlmax[ir2];
        int kkr1=(lmax1+1)*(lmax1+1);
        int kkr2=(lmax2+1)*(lmax2+1);
        int kkr1_ns = kkr1 * lsms.n_spin_cant;
        int kkr2_ns = kkr2 * lsms.n_spin_cant;
        Real rij[3];
        rij[0]=atom.LIZPos(0,ir1)-atom.LIZPos(0,ir2);
        rij[1]=atom.LIZPos(1,ir1)-atom.LIZPos(1,ir2);
        rij[2]=atom.LIZPos(2,ir1)-atom.LIZPos(2,ir2);
        buildBGijCPU(lsms, atom, ir1, ir2, rij, energy, prel, iOffset, jOffset, bgij);

#ifdef COMPARE_ORIGINAL
        int kkri=(lmax1+1)*(lmax1+1);
        int kkrj=(lmax2+1)*(lmax2+1);
        int lmax=lsms.maxlmax;
        int kkrsz=(lmax+1)*(lmax+1);
        bool exitCompare = false;
        Matrix<Complex> tmat_n(lsms.n_spin_cant*atom.kkrsz, lsms.n_spin_cant*atom.kkrsz);
        int im=0;
        if(lsms.n_spin_pola == lsms.n_spin_cant) // non polarized or spin canted
        {
          for(int js=0; js<lsms.n_spin_cant; js++)
          {
            int jsm = kkrsz*kkrsz_ns*js;
            for(int j=0; j<kkr1; j++)
            {
              for(int is=0; is<lsms.n_spin_cant; is++)
              {
                int jm=jsm+kkrsz_ns*j+kkrsz*is;
                int one=1;
                BLAS::zcopy_(&kkr1,&local.tmatStore(iie*local.blkSizeTmatStore+jm,atom.LIZStoreIdx[ir1]),&one,&tmat_n[im],&one);
                im+=kkr1;
              }
            }
          }
        } else { // spin polarized colinear version for ispin
          int ispin=0;
          printf("warning: cant't test building kkrMatrix for collinear spin polarized yet!\n");
          exit(1);
          int jsm = kkrsz*kkrsz*ispin; // copy spin up or down?
          for(int j=0; j<kkr1; j++)
          {
            int jm=jsm+kkrsz_ns*j;
            int one=1;
            BLAS::zcopy_(&kkr1,&local.tmatStore(iie*local.blkSizeTmatStore+jm,atom.LIZStoreIdx[ir1]),&one,&tmat_n[im],&one);
            im+=kkr1;
          }
        }

        for(int i=0; i<2*kkri; i++)
          for(int j=0; j<2*kkrj; j++)
          {
            int jm = i + j*2*kkri;
            if(tmat_n(i, j) != local.tmatStore(iie*local.blkSizeTmatStore+jm,atom.LIZStoreIdx[ir1]))
            {
              printf("buildKKRMatrixLMaxDifferentCPU: tmat_n(%d, %d) [%g + %gi] != tmatStore(%d + %d, %d) [%g + %gi]\n",
                     i, j, tmat_n(i, j).real(), tmat_n(i, j).imag(),
                     iie*local.blkSizeTmatStore, jm,atom.LIZStoreIdx[ir1],
                     local.tmatStore(iie*local.blkSizeTmatStore+jm,atom.LIZStoreIdx[ir1]).real(),
                     local.tmatStore(iie*local.blkSizeTmatStore+jm,atom.LIZStoreIdx[ir1]).imag());
              exitCompare = true;
            }
          }

        if(exitCompare) exit(1);
        
        Matrix<Complex> bgijTest(kkr1_ns, kkr2_ns);
        for(int i=0; i<kkr1_ns; i++)
          for(int j=0; j<kkr2_ns; j++)
            bgijTest(i,j) = bgij(iOffset + i, jOffset + j);

        /*
        BLAS::zgemm_("n", "n", &kkr1_ns, &kkr2_ns, &kkr1_ns, &cmone,
                     &tmat_n(0, 0), &kkr1_ns,
                     &bgijTest(0, 0), &kkr1_ns, &czero,
                     &m(iOffset, jOffset), &nrmat_ns);
        */
#endif
        if(lsms.n_spin_pola == lsms.n_spin_cant) // non polarized or spin canted
        {
          BLAS::zgemm_("n", "n", &kkr1_ns, &kkr2_ns, &kkr1_ns, &cmone,
                       &local.tmatStore(iie*local.blkSizeTmatStore, atom.LIZStoreIdx[ir1]), &kkrsz_ns,
                       &bgij(iOffset, jOffset), &nrmat_ns, &czero,
                       &m(iOffset, jOffset), &nrmat_ns);
        } else {  // spin polarized, collinear
          int lmax=lsms.maxlmax;
          int kkrsz=(lmax+1)*(lmax+1);
          int spinOffset = kkrsz * kkrsz * ispin;
          BLAS::zgemm_("n", "n", &kkr1_ns, &kkr2_ns, &kkr1_ns, &cmone,
                       &local.tmatStore(iie*local.blkSizeTmatStore + spinOffset, atom.LIZStoreIdx[ir1]), &kkrsz_ns,
                       &bgij(iOffset, jOffset), &nrmat_ns, &czero,
                       &m(iOffset, jOffset), &nrmat_ns);
        }
       
      }
    }
  }
#ifdef COMPARE_ORIGINAL
  bool exitCompare = false;
  int ispin = 0;
  Matrix<Complex> mTest(nrmat_ns, nrmat_ns);
  buildKKRMatrix(lsms, local, atom, ispin, energy, prel, iie, mTest);
  for(int i=0; i<nrmat_ns; i++)
    for(int j=0; j<nrmat_ns; j++)
    {
      if(m(i, j) != mTest(i,j))
      {
        printf("buildKKRMatrixLMaxDifferentCPU: m(%d, %d) [%g + %gi] != mTest(%d, %d) [%g + %gi]\n",
               i, j, m(i, j).real(), m(i, j).imag(),
               i, j, mTest(i, j).real(), mTest(i, j).imag());
        exitCompare = true;
      }
    }
  
  if(exitCompare) exit(1);
#endif
}

void buildKKRMatrixCPU(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int ispin, int iie, Complex energy, Complex prel,
                                    Matrix<Complex> &m)
{
  // decide between identical lmax and different lmax:
  
  bool lmaxIdentical = true;

  if(atom.LIZlmax[0] != lsms.maxlmax)
  {
    lmaxIdentical = false;
    printf("atom.LIZlmax[0] (=%d) != lsms.maxlmax (=%d)\n",atom.LIZlmax[0], lsms.maxlmax);
  }
  for(int ir = 0; ir < atom.numLIZ; ir++)
  {
    if(atom.LIZlmax[ir] != atom.LIZlmax[0])
      lmaxIdentical = false;
  }
  
  if(lmaxIdentical)
  {
    // printf("lmax identical in buildKKRMatrix\n");
    buildKKRMatrixLMaxIdenticalCPU(lsms, local, atom, ispin, iie, energy, prel, m);
  } else {
    // printf("lmax not identical in buildKKRMatrix\n");
    buildKKRMatrixLMaxDifferentCPU(lsms, local, atom, ispin, iie, energy, prel, m);
  }
}
