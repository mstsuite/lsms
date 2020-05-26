/* -*- mode: C++; c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#include "buildKKRMatrix.hpp"

#include <stdio.h>

#include "Complex.hpp"
#include "Matrix.hpp"
#include <vector>

#include "Accelerator/DeviceStorage.hpp"
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

// we might want to distinguish between systems where all lmax (and consequently kkrsz_ns) are the same
// and systems with potential different lmax on different atoms and l steps

// Fortran layout for matrix
#define IDX(i, j, lDim) (((j)*(lDim))+(i))

__device__
inline void calculateHankelCuda(cuDoubleComplex prel, double r, int lend, cuDoubleComplex *hfn)
{
  if(threadIdx.x == 0)
  {
    const cuDoubleComplex sqrtm1(0.0, 1.0);
    cuDoubleComplex z=prel*r;
    hfn[0]=-sqrtm1;
    hfn[1]=-1.0-sqrtm1/z;
    for(int l=1; l<lend; l++)
    {
      hfn[l+1]=(2.0*l+1)*hfn[l]/z - hfn[l-1];
    }
  /*
c             l+1
c     hfn = -i   *h (k*R  )*sqrt(E)
c                  l    ij
*/
    z=std::exp(sqrtm1*z)/rmag;
    for(int l=0; l<=lend;l++)
    {
      hfn[l]=-hfn[l]*z*ilp1[l]; 
    }
  }
//  __syncthreads();
}

__device__
inline void calculateSinCosPowersCuda(Real *rij, int lend, Real *sinmp, Real *cosmp)
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

__device__ __inline__ int plmIdxDev(int l, int m)
{ return l*(l+1)/2+m; }

__device__
void associatedLegendreFunctionNormalizedCuda(Real x, int lmax, Real *Plm)
{
  const Real pi = std::acos(-R(1));
  // y = \sqrt{1-x^2}
  Real y = std::sqrt(1.0-x*x);
  // initialize the first entry
  // Plm[0]=std::sqrt(R(1)/(R(2)*pi));
  Plm[0]=std::sqrt(1.0/(4.0*pi));

  if(lmax<1) return;

  for(int m=1; m<=lmax; m++)
  {
    // \bar{P}_{mm} = - \sqrt{\frac{2m+1}{2m}} y \bar{P}_{m-1, m-1}
    Plm[plmIdxDev(m,m)] = - std::sqrt(Real(2*m+1)/Real(2*m)) * y * Plm[plmIdxDev(m-1,m-1)];
    // \bar{P}_{mm-1} = \sqrt{2 m + 1} x \bar{P}_{m-1, m-1}
    Plm[plmIdxDev(m,m-1)] = std::sqrt(Real(2*m+1)) * x * Plm[plmIdxDev(m-1,m-1)]; 
  }

  for(int m=0; m<lmax; m++)
  {
    for(int l=m+2; l<=lmax; l++)
    {
      // \bar{P}_{lm} = a_{lm} (x \bar{P}_{l-1. m} - b_{lm} \bar{P}_{l-2, m})
      // a_{lm} = \sqrt{\frac{(4 l^2 - 1)(l^2 - m^2)}}
      // b_{lm} = \sqrt{\frac{(l -1)^2 - m^2}{4 (l-1)^2 -1}}
      R a_lm = std::sqrt(Real(4*l*l-1)/Real(l*l - m*m));
      R b_lm = std::sqrt(Real((l-1)*(l-1) - m*m)/Real(4*(l-1)*(l-1)-1));
      Plm[plmIdxDev(l,m)] = a_lm * (x * Plm[plmIdxDev(l-1,m)] - b_lm * Plm[plmIdxDev(l-2,m)]);
    }
  }
}

size_t sharedMemoryBGijCuda(LSMSSystemParameters &lsms, size_t *hfnOffset, size_t *sinmpOffset, size_t *cosmpOffset,
                            size_t *plmOffset, size_t *dlmOffset)
{
  size = 0;

  *hfnOffset = size;
  size += sizeof(cuDoubleComplex) * (2*lsms.maxlmax + 1);

  *sinmpOffset = size;
  size += sizeof(double) * (2*lsms.maxlmax + 1);

  *cosmpOffset = size;
  size += sizeof(double) * (2*lsms.maxlmax + 1);

  *plmOffset = size;
  size += sizeof(double) * (lsms.angularMomentumIndices.ndlm);

  *dlmOffset = size;
  size += sizeof(cuDoubleComplex) * (lsms.angularMomentumIndices.ndlj);
  
  return size;
}

__global__
void buildBGijCudaKernel(Real *LIZPos, int *LIZlmax, int *lofk, int *mofk,
                         size_t hfnOffset, size_t sinmpOffset, size_t cosmpOffset, size_t plmOffset, size_t dlmOffset,
                         cuDoubleComplex energy, cuDoubleComplex prel, int *offsets, cuDoubleComplex *devBgij)
//  void buildBGijCPU(LSMSSystemParameters &lsms, AtomData &atom, int ir1, int ir2, Real *rij,
//                  Complex energy, Complex prel, int iOffset, int jOffset, Matrix<Complex> &bgij)
{
  int ir1 = blockIdx.x;
  int ir2 = blockIdx.y;
  extern char __shared__ sharedMemory[]; 

  if(ir1 != ir2)
  {
    int iOffset = offsets[ir1];
    // int iOffset = ir1 * kkrsz_ns;
    int jOffset = offsets[ir2];
    // int jOffset = ir2 * kkrsz_ns;

    Real rij[3];
    rij[0] = LIZPos[3*ir1 + 0] - LIZPos[3*ir2 + 0];
    rij[1] = LIZPos[3*ir1 + 1] - LIZPos[3*ir2 + 1];
    rij[2] = LIZPos[3*ir1 + 2] - LIZPos[3*ir2 + 2];
  
    // Complex hfn[2*lsms.maxlmax + 1];
    cuDobleComplex *hfn = (cuDoubleComplex *) (sharedMemory + hfnOffset);
    // Real sinmp[2*lsms.maxlmax + 1];
    Real *sinmp = (Real *) (sharedMemory + sinmpOffset);
    // Real cosmp[2*lsms.maxlmax + 1];
    Real *cosmp = (Real *) (sharedMemory + cosmpOffset);
    // Real plm[lsms.angularMomentumIndices.ndlm];
    Real *plm = (Real *) (sharedMemory + plmOffset);
    // Complex dlm[lsms.angularMomentumIndices.ndlj];
    cuDoubleComplex *dlm = (cuDoubleComplex + dlmOffset);
    Real r = std::sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
    int lmax1 = LIZlmax[ir1];
    int lmax2 = LIZlmax[ir2];
    int kkri=(lmax1+1)*(lmax1+1);
    int kkrj=(lmax2+1)*(lmax2+1);
    int lend = lmax1 + lmax2;

    Real pi4=4.0*2.0*std::asin(1.0);
    Real cosTheta = rij[2]/r;
  
    if(threadIdx.x == 0)
    {
      calculateHankelCuda(prel, r, lend, hfn);

      associatedLegendreFunctionNormalizedCuda(cosTheta, lend, plm);
      // for associatedLegendreFunctionNormalized all clm[i] == 1.0
      // for(int j=0;j<ndlm_local;j++)
      //   plm[j]=clm[j]*plm[j];
  
      //     calculate cos(phi) and sin(phi) .................................
      // needs to be serial
      calculateSinCosPowersCuda(rij, lend, sinmp, cosmp);
    }
    __syncthreads();
  
    // can be parallel
    int j=0;
    for(int l = threadIdx.x; l<=lend; l += blockDim.x)
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

    // for(int i=0; i<kkri; i++)
    //   for(int j=0; j<kkrj; j++)
    // for(int ij=0; ij < kkri*kkrj; ij++)
    for(int ij=threadIdx.x; ij < kkri*kkrj; ij += blockDim.x)
    {
      int lm2 = ij % kkri;
      int lm1 = ij / kkri;
      devBgi[IDX(iOffset + lm2, jOffset + lm1, nrmat_ns)] = 0.0;
      // bgij(iOffset + lm2, jOffset + lm1) = 0.0;
      // }
    
//     loop over l1,m1............................................
//    for(int lm1=0; lm1<kkrj; lm1++)
//    {
      int l1=lofk[lm1];
      int m1=mofk[lm1];
    
//        loop over l2,m2..............................................
//      for(int lm2=0; lm2<kkri; lm2++)
//      {
      int l2=lofk[lm2];
      int m2=mofk[lm2];
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
        devBgij[IDX(iOffset + lm2, jOffset + lm1, nrmat_ns)] += GauntCoeficients::cgnt(l3/2,lm1,lm2)*dlm[j];
      }
      // gij[lm2+lm1*kkri]=pi4*illp(lm2,lm1)*gij[lm2+lm1*kkri];
      devBgij[IDX(iOffset + lm2, jOffset + lm1, nrmat_ns)] *= pi4 * IFactors::illp(lm2,lm1);
    }

  // might do this as a seperate kernel
    __syncthread();
    setBGijCuda(lsms, atom, ir1, ir2, iOffset, jOffset, bgij);
  }
}

void buildKKRMatrixLMaxIdenticalCuda(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom, int iie,
                                     Complex energy, Complex prel,
                                     Complex *tMatrix, Complex *devM)
{
  cublasHandle_t cublasHandle = DeviceStorage::getCublasHandle();
  int nrmat_ns = lsms.n_spin_cant*atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block

  Complex cmone = Complex(-1.0,0.0);
  Complex czero=0.0;

  Complex *devBgij = d.getDevBGij();
  // Matrix<Complex> bgijSmall(kkrsz_ns, kkrsz_ns);

  unitMatrixCuda<Complex>(devM, nrmat_ns, nrmat_ns);
  zeroMatrixCuda(devBgij, nrmat_ns, nrmat_ns);

// calculate Bgij
// reuse ipvt for offsets
  int *devOffsets = d.getDevIpvt();
  
  std::vector<int> offsets(atom.numLIZ);
  for(int ir = 0; ir < atom.numLIZ; ir++)
    offsets[ir] = ir * kkrsz_ns;

  cudaMemcpy(devOffsets, &offsets[0], atom.numLIZ*sizeof(int), cudaMemcpyHostToDevice);

  size_t hfnOffset, sinmpOffset, cosmpOffset, plmOffset, dlmOffset;
  size_t smSize = sharedMemoryBGijCuda(lsms, &hfnOffset, &sinmpOffset, &cosmpOffset,
                                       &plmOffset, &dlmOffset);
  int treads = 256;
  dim3 blocks = dim3(atom.numLIZ,atom.numLIZ,1);
  buildBGijCudaKernel<<<blocks,treads,smSize>>>(devAtom.LIZPos, devAtom.LIZlmax,
                                                devConstants.lofk, devConstants.mofk,
                                                hfnOffset, sinmpOffset, cosmpOffset, plmOffset, dlmOffset,
                                                energy, prel, devOffsets, devBgij);
  
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
        
        buildBGijCuda(lsms, atom, ir1, ir2, rij, energy, prel, iOffset, jOffset, bgij);
        // buildBGijCPU(lsms, atom, ir1, ir2, rij, energy, prel, 0, 0, bgijSmall);
             
        BLAS::zgemm_("n", "n", &kkr1_ns, &kkr2_ns, &kkr1_ns, &cmone,
                     &local.tmatStore(iie*local.blkSizeTmatStore, atom.LIZStoreIdx[ir1]), &kkr1_ns,
                     // &tmat_n(0, 0), &kkr1_ns,
                     &bgij(iOffset, jOffset), &nrmat_ns, &czero,
                     // &bgijSmall(0, 0), &kkrsz_ns, &czero,
                     &m(iOffset, jOffset), &nrmat_ns);
        
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
}

void buildKKRMatrixLMaxIdenticalCuda(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom, int iie,
                        Complex *tMatrix, Complex *devM)
{
  cublasHandle_t cublasHandle = DeviceStorage::getCublasHandle();
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
        
        buildBGijCuda(lsms, atom, ir1, ir2, rij, energy, prel, iOffset, jOffset, bgij);
        // buildBGijCPU(lsms, atom, ir1, ir2, rij, energy, prel, 0, 0, bgijSmall);
             
        BLAS::zgemm_("n", "n", &kkr1_ns, &kkr2_ns, &kkr1_ns, &cmone,
                     &local.tmatStore(iie*local.blkSizeTmatStore, atom.LIZStoreIdx[ir1]), &kkr1_ns,
                     // &tmat_n(0, 0), &kkr1_ns,
                     &bgij(iOffset, jOffset), &nrmat_ns, &czero,
                     // &bgijSmall(0, 0), &kkrsz_ns, &czero,
                     &m(iOffset, jOffset), &nrmat_ns);
        
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
}

void buildKKRMatrixCuda(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie, Complex energy, Complex prel,
                       Complex *devM)
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
    buildKKRMatrixLMaxIdenticalCuda(lsms, local, atom, iie, energy, prel, m);
  } else {
    // printf("lmax not identical in buildKKRMatrix\n");
     buildKKRMatrixLMaxDifferentCuda(lsms, local, atom, iie, energy, prel, m);
  }
}
