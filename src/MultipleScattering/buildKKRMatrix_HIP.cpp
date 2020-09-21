/* -*- mode: C++; c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#include "buildKKRMatrix.hpp"

#include <stdio.h>

#include "Complex.hpp"
#include "Matrix.hpp"
#include <vector>

#include "Accelerator/DeviceStorage.hpp"
#include <hip/hip_runtime.h>
#include <hip/hip_complex.h>
#include <hipblas.h>
#include <rocsolver.h>

#include "linearSolvers.hpp"

// we might want to distinguish between systems where all lmax (and consequently kkrsz_ns) are the same
// and systems with potential different lmax on different atoms and l steps

// Fortran layout for matrix
#define IDX(i, j, lDim) (((j)*(lDim))+(i))
#define IDX3(i, j, k, lDim, mDim) (((k)*(lDim)*(mDim)) + ((j)*(lDim)) + (i))

__device__ __inline__
deviceDoubleComplex complexExp(deviceDoubleComplex z)
{
  double mExp=exp(hipCreal(z));
  return make_hipDoubleComplex(mExp*cos(hipCimag(z)), mExp*sin(hipCimag(z)));
}

__device__
inline void calculateHankelHip(deviceDoubleComplex prel, double r, int lend, deviceDoubleComplex *ilp1, deviceDoubleComplex *hfn)
{
  if(hipThreadIdx_x == 0)
  {
    const deviceDoubleComplex sqrtm1 = make_hipDoubleComplex(0.0, 1.0);
    deviceDoubleComplex z=prel*make_hipDoubleComplex(r,0.0);
    hfn[0]=make_hipDoubleComplex(0.0, -1.0); //-sqrtm1;
    hfn[1]=-1.0-sqrtm1/z;
    for(int l=1; l<lend; l++)
    {
      hfn[l+1]=(2.0*l+1)*hfn[l]/z - hfn[l-1];
    }

//             l+1
//     hfn = -i   *h (k*R  )*sqrt(E)
//                  l    ij

    z=complexExp(sqrtm1*z)/r;
    for(int l=0; l<=lend;l++)
    {
      hfn[l]=-hfn[l]*z*ilp1[l]; 
    }
  }
//  __syncthreads();
}


__device__
inline void calculateSinCosPowersHip(Real *rij, int lend, Real *sinmp, Real *cosmp)
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
void associatedLegendreFunctionNormalizedHip(Real x, int lmax, Real *Plm)
{
  const Real pi = std::acos(-1.0);
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
      Real a_lm = std::sqrt(Real(4*l*l-1)/Real(l*l - m*m));
      Real b_lm = std::sqrt(Real((l-1)*(l-1) - m*m)/Real(4*(l-1)*(l-1)-1));
      Plm[plmIdxDev(l,m)] = a_lm * (x * Plm[plmIdxDev(l-1,m)] - b_lm * Plm[plmIdxDev(l-2,m)]);
    }
  }
}


size_t sharedMemoryBGijHip(LSMSSystemParameters &lsms, size_t *hfnOffset, size_t *sinmpOffset, size_t *cosmpOffset,
                            size_t *plmOffset, size_t *dlmOffset)
{
  size_t size = 0;

  *hfnOffset = size;
  size += sizeof(deviceDoubleComplex) * (2*lsms.maxlmax + 1);

  *sinmpOffset = size;
  size += sizeof(double) * (2*lsms.maxlmax + 1);

  *cosmpOffset = size;
  size += sizeof(double) * (2*lsms.maxlmax + 1);

  *plmOffset = size;
  size += sizeof(double) * (lsms.angularMomentumIndices.ndlm);

  *dlmOffset = size;
  size += sizeof(deviceDoubleComplex) * (lsms.angularMomentumIndices.ndlj);
  
  return size;
}

__global__
void setBGijHip(bool fullRelativity, int n_spin_cant, int *LIZlmax,
                 int *offsets, size_t nrmat_ns, deviceDoubleComplex *devBgij)
{
  if(n_spin_cant == 1) return;

  int ir1 = hipBlockIdx_x;
  int ir2 = hipBlockIdx_y;
  int iOffset = offsets[ir1];
  int jOffset = offsets[ir2];

  int kkri=(LIZlmax[ir1]+1)*(LIZlmax[ir1]+1);
  int kkrj=(LIZlmax[ir2]+1)*(LIZlmax[ir2]+1);

  if(!fullRelativity) //(lsms.relativity != full)
  {
    for(int ij=hipThreadIdx_x; ij < kkri*kkrj; ij += hipBlockDim_x)
    {
      int i = ij % kkri;
      int j = ij / kkri;

      devBgij[IDX(iOffset + kkri + i, jOffset        + j, nrmat_ns)] = make_hipDoubleComplex(0.0, 0.0); // bgij(iOffset + i, jOffset + j);
      devBgij[IDX(iOffset        + i, jOffset + kkrj + j, nrmat_ns)] = make_hipDoubleComplex(0.0, 0.0); // bgij(iOffset + i, jOqffset + j);
      devBgij[IDX(iOffset + kkri + i, jOffset + kkrj + j, nrmat_ns)] = devBgij[IDX(iOffset + i, jOffset + j, nrmat_ns)];
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
    // exit(1);
  }
}


__global__
void buildGijHipKernel(Real *LIZPos, int *LIZlmax, int *lofk, int *mofk, deviceDoubleComplex *ilp1, deviceDoubleComplex *illp, Real *cgnt,
                        int ndlj_illp, int lmaxp1_cgnt, int ndlj_cgnt,
                         size_t hfnOffset, size_t sinmpOffset, size_t cosmpOffset, size_t plmOffset, size_t dlmOffset,
                         deviceDoubleComplex energy, deviceDoubleComplex prel, int *offsets, size_t nrmat_ns, deviceDoubleComplex *devBgij)
//  void buildBGijCPU(LSMSSystemParameters &lsms, AtomData &atom, int ir1, int ir2, Real *rij,
//                  Complex energy, Complex prel, int iOffset, int jOffset, Matrix<Complex> &bgij)
{
  int ir1 = hipBlockIdx_x;
  int ir2 = hipBlockIdx_y;
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
    deviceDoubleComplex *hfn = (deviceDoubleComplex *) (sharedMemory + hfnOffset);
    // Real sinmp[2*lsms.maxlmax + 1];
    Real *sinmp = (Real *) (sharedMemory + sinmpOffset);
    // Real cosmp[2*lsms.maxlmax + 1];
    Real *cosmp = (Real *) (sharedMemory + cosmpOffset);
    // Real plm[lsms.angularMomentumIndices.ndlm];
    Real *plm = (Real *) (sharedMemory + plmOffset);
    // Complex dlm[lsms.angularMomentumIndices.ndlj];
    deviceDoubleComplex *dlm = (deviceDoubleComplex *) (sharedMemory + dlmOffset);
    Real r = std::sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
    int lmax1 = LIZlmax[ir1];
    int lmax2 = LIZlmax[ir2];
    int kkri=(lmax1+1)*(lmax1+1);
    int kkrj=(lmax2+1)*(lmax2+1);
    int lend = lmax1 + lmax2;

    Real pi4=4.0*2.0*std::asin(1.0);
    Real cosTheta = rij[2]/r;
  
    if(hipThreadIdx_x == 0)
    {
      calculateHankelHip(prel, r, lend, ilp1, hfn);

      associatedLegendreFunctionNormalizedHip(cosTheta, lend, plm);
      // for associatedLegendreFunctionNormalized all clm[i] == 1.0
      // for(int j=0;j<ndlm_local;j++)
      //   plm[j]=clm[j]*plm[j];
  
      //     calculate cos(phi) and sin(phi) .................................
      // needs to be serial
      calculateSinCosPowersHip(rij, lend, sinmp, cosmp);
    }
    __syncthreads();
  
    // can be parallel
    int j=0;
    for(int l = hipThreadIdx_x; l<=lend; l += hipBlockDim_x)
    {
      int ll = l*(l+1);
      j = ll;
      ll = ll/2;
      Real m1m = 1.0;
      dlm[j] = hfn[l]*plm[ll];
      for(int m=1; m<=l; m++)
      {
        m1m = -m1m;
        deviceDoubleComplex fac = plm[ll+m] * make_hipDoubleComplex(cosmp[m],sinmp[m]);
        dlm[j-m] = hfn[l] * m1m * fac;
        dlm[j+m] = hfn[l] * hipConj(fac);
      }
    }
//     ================================================================
//     calculate g(R_ij)...............................................
  // for(int i=0; i<kkri*kkrj; i++) gij[i]=0.0;

    // for(int i=0; i<kkri; i++)
    //   for(int j=0; j<kkrj; j++)
    // for(int ij=0; ij < kkri*kkrj; ij++)
    for(int ij = hipThreadIdx_x; ij < kkri*kkrj; ij += hipBlockDim_x)
    {
      int lm2 = ij % kkri;
      int lm1 = ij / kkri;
      devBgij[IDX(iOffset + lm2, jOffset + lm1, nrmat_ns)] = make_hipDoubleComplex(0.0, 0.0);
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
        
//          ==========================================================
//          l2-l1
//          illp(lm2,lm1) = i
//
//          perform sum over l3 with gaunt # ......................
//          ==========================================================
        
      int m3=m2-m1;
      int llow=max(abs(m3), abs(l1-l2));
      if(hipCabs(prel)==0.0) llow=l1+l2;
      for(int l3=l1+l2; l3>=llow; l3-=2)
      {
        int j=l3*(l3+1)+m3;
        // gij[lm2+lm1*kkri] = gij[lm2+lm1*kkri]+cgnt(l3/2,lm1,lm2)*dlm[j];
        devBgij[IDX(iOffset + lm2, jOffset + lm1, nrmat_ns)] =  devBgij[IDX(iOffset + lm2, jOffset + lm1, nrmat_ns)]
          + cgnt[IDX3(l3/2,lm1,lm2,lmaxp1_cgnt,ndlj_cgnt)]*dlm[j];
      }
      // gij[lm2+lm1*kkri]=pi4*illp(lm2,lm1)*gij[lm2+lm1*kkri];
      devBgij[IDX(iOffset + lm2, jOffset + lm1, nrmat_ns)] = devBgij[IDX(iOffset + lm2, jOffset + lm1, nrmat_ns)]
        * pi4 * illp[IDX(lm2, lm1, ndlj_illp)];
    }

  // might do this as a seperate kernel
  //  __syncthreads();
  //  setBGijCuda(fullRelativity, n_spin_cant, LIZlmax, ir1, ir2, iOffset, jOffset, nrmat_ns, devBgij);
  }
}


__device__
void buildTmatNHip(int ispin, int n_spin_pola, int n_spin_cant, int iie, int blkSizeTmatStore, int tmatStoreLDim,
                    int kkr1, int kkr2, int lizStoreIdx,
                    deviceDoubleComplex *devTmatStore, int kkrsz_ns, deviceDoubleComplex *tmat_n)
{
// Matrix<Complex> tmat_n(lsms.n_spin_cant*atom.kkrsz, lsms.n_spin_cant*atom.kkrsz);
  if(hipThreadIdx_x == 0)
  {
    int im=0;
    if(n_spin_pola == n_spin_cant) // non polarized or spin canted
    {
      int kkrsz = kkrsz_ns/n_spin_cant;
      for(int js=0; js<n_spin_cant; js++)
      {
        int jsm = kkrsz*kkrsz_ns*js;
        for(int j=0; j<kkr1; j++)
        {
          for(int is=0; is<n_spin_cant; is++)
          {
            int jm=jsm+kkrsz_ns*j+kkrsz*is;
//                int one=1;
//                BLAS::zcopy_(&kkr1,&local.tmatStore(iie*local.blkSizeTmatStore+jm,atom.LIZStoreIdx[ir1]),&one,&tmat_n[im],&one);
            for(int i=0; i<kkr1; i++)
            {
              tmat_n[im+i] = devTmatStore[IDX(iie*blkSizeTmatStore+jm+i, lizStoreIdx, tmatStoreLDim)];
            }
            im+=kkr1;
          }
        }
      }
    } else { // spin polarized colinear version for ispin
      int kkrsz = kkrsz_ns/n_spin_cant; 
      // int ispin=0;
      printf("warning: cant't test building kkrMatrix for collinear spin polarized yet!\n");
      // exit(1);
      int jsm = kkrsz*kkrsz*ispin; // copy spin up or down?
      for(int j=0; j<kkr1; j++)
      {
        int jm=jsm+kkrsz_ns*j;
//           int one=1;
//           BLAS::zcopy_(&kkr1,&local.tmatStore(iie*local.blkSizeTmatStore+jm,atom.LIZStoreIdx[ir1]),&one,&tmat_n[im],&one);
        for(int i=0; i<kkr1; i++)
        {
          tmat_n[im+i] = devTmatStore[IDX(iie*blkSizeTmatStore+jm+i, lizStoreIdx, tmatStoreLDim)];
        }
        im+=kkr1;
      }
    }
  }
  __syncthreads();
}

__global__
void buildKKRMatrixMultiplyKernelHip(int *LIZlmax, int *LIZStoreIdx, int *offsets, int kkrsz_ns,
                                      int ispin, int n_spin_pola, int n_spin_cant, int iie, int blkSizeTmatStore, int tmatStoreLDim,
                                      deviceDoubleComplex *devTmatStore, int nrmat_ns, deviceDoubleComplex *devBgij, deviceDoubleComplex *devM)
{
  int ir1 = blockIdx.x;
  int ir2 = blockIdx.y;
  deviceDoubleComplex __shared__ *tmat_n;
  int iOffset = offsets[ir1];
  int jOffset = offsets[ir2];

  if(ir1 != ir2)
  {
    int lmax1 = LIZlmax[ir1];
    int lmax2 = LIZlmax[ir2];
    int kkr1=(lmax1+1)*(lmax1+1);
    int kkr2=(lmax2+1)*(lmax2+1);
    int kkr1_ns = kkr1 * n_spin_cant;
    int kkr2_ns = kkr2 * n_spin_cant;
        
//        BLAS::zgemm_("n", "n", &kkr1_ns, &kkr2_ns, &kkr1_ns, &cmone,
//                     &local.tmatStore(iie*local.blkSizeTmatStore, devAtom.LIZStoreIdx[ir1]), &kkr1_ns,
//                     // &tmat_n(0, 0), &kkr1_ns,
//                     &bgij(iOffset, jOffset), &nrmat_ns, &czero,
//                     // &bgijSmall(0, 0), &kkrsz_ns, &czero,
//                     &m(iOffset, jOffset), &nrmat_ns);
        
//    for(int i=0; i<kkr1_ns; i++)
//      for(int j=0; j<kkr2_ns; j++)
    buildTmatNHip(ispin, n_spin_pola, n_spin_cant, iie, blkSizeTmatStore, tmatStoreLDim,
                   kkr1, kkr2, LIZStoreIdx[ir1], devTmatStore, kkrsz_ns, tmat_n);
    
    for(int ij=hipThreadIdx_x; ij < kkr1_ns*kkr2_ns; ij += hipBlockDim_x)
    {
      int i = ij % kkr1_ns;
      int j = ij / kkr1_ns;

      devM[IDX(iOffset + i, jOffset + j, nrmat_ns)] = make_hipDoubleComplex(0.0,0.0);
      for(int k=0; k<kkr1_ns ; k++)
        devM[IDX(iOffset + i, jOffset + j, nrmat_ns)] = devM[IDX(iOffset + i, jOffset + j, nrmat_ns)] -
          tmat_n[IDX(i,j,kkrsz_ns)] * // tmat_n(i, k) * // local.tmatStore(iie*local.blkSizeTmatStore + , atom.LIZStoreIdx[ir1]) *
          devBgij[IDX(iOffset + k, jOffset + j, nrmat_ns)];
    }
    
  }
}


/*
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
        
//        for(int i=0; i<kkr1_ns; i++)
//          for(int j=0; j<kkr2_ns; j++)
//          {
//            m(iOffset + i, jOffset + j) = 0.0;
//            for(int k=0; k<kkr1_ns ; k++)
//              m(iOffset + i, jOffset + j) -= tmat_n(i, k) * // local.tmatStore(iie*local.blkSizeTmatStore + , atom.LIZStoreIdx[ir1]) *
//                // bgij(iOffset + k, jOffset + j);
//                bgijSmall(k, j);
//          }
        
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
        
        
//        for(int i=0; i<kkr1_ns; i++)
//          for(int j=0; j<kkr2_ns; j++)
//          {
//            m(iOffset + i, jOffset + j) = 0.0;
//            for(int k=0; k<kkr1_ns ; k++)
//              m(iOffset + i, jOffset + j) -= tmat_n(i, k) * // local.tmatStore(iie*local.blkSizeTmatStore + , atom.LIZStoreIdx[ir1]) *
//                // bgij(iOffset + k, jOffset + j);
//                bgijSmall(k, j);
//          }
        
      }
    }
  }
}

*/

void buildKKRMatrixLMaxDifferentHip(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom,
                                     DeviceStorage &d, DeviceAtom &devAtom, int ispin,
                                     int iie, Complex energy, Complex prel, Complex *devM)
{
  hipblasHandle_t hipblasHandle = DeviceStorage::getHipBlasHandle();
  int nrmat_ns = lsms.n_spin_cant*atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block
  bool fullRelativity = false;
  if(lsms.relativity == full) fullRelativity = true;

  // Complex cmone = Complex(-1.0,0.0);
  // Complex czero=0.0;

  Complex *devBgij = d.getDevBGij();
  // Matrix<Complex> bgijSmall(kkrsz_ns, kkrsz_ns);

  deviceDoubleComplex hipEnergy = make_hipDoubleComplex(energy.real(), energy.imag());
  deviceDoubleComplex hipPrel = make_hipDoubleComplex(prel.real(), prel.imag());
  
  unitMatrixHip<Complex>(devM, nrmat_ns, nrmat_ns);
  zeroMatrixHip(devBgij, nrmat_ns, nrmat_ns);

// calculate Bgij
// reuse ipvt for offsets
  int *devOffsets = d.getDevIpvt();
  
  std::vector<int> offsets(devAtom.numLIZ);
  offsets[0] = 0;
  for(int ir = 1; ir < atom.numLIZ; ir++)
    offsets[ir] = offsets[ir-1] + lsms.n_spin_cant * (atom.LIZlmax[ir-1]+1)*(atom.LIZlmax[ir-1]+1);

  deviceMemcpy(devOffsets, &offsets[0], atom.numLIZ*sizeof(int), deviceMemcpyHostToDevice);

  size_t hfnOffset, sinmpOffset, cosmpOffset, plmOffset, dlmOffset;
  size_t smSize = sharedMemoryBGijHip(lsms, &hfnOffset, &sinmpOffset, &cosmpOffset,
                                       &plmOffset, &dlmOffset);
  int threads = 256;
  dim3 blocks = dim3(devAtom.numLIZ, devAtom.numLIZ,1);
  buildGijHipKernel<<<blocks,threads,smSize>>>(devAtom.LIZPos, devAtom.LIZlmax,
                                                DeviceConstants::lofk, DeviceConstants::mofk, DeviceConstants::ilp1, DeviceConstants::illp, DeviceConstants::cgnt,
                                                DeviceConstants::ndlj_illp, DeviceConstants::lmaxp1_cgnt, DeviceConstants::ndlj_cgnt,
                                                hfnOffset, sinmpOffset, cosmpOffset, plmOffset, dlmOffset,
                                                hipEnergy, hipPrel,
                                                devOffsets, nrmat_ns, (deviceDoubleComplex *)devBgij);

  setBGijHip<<<blocks, threads>>>(fullRelativity, lsms.n_spin_cant, devAtom.LIZlmax,
                                   devOffsets, nrmat_ns, (deviceDoubleComplex *)devBgij);


  smSize = kkrsz_ns*kkrsz_ns*sizeof(deviceDoubleComplex);
  buildKKRMatrixMultiplyKernelHip<<<blocks, threads, smSize>>>(devAtom.LIZlmax, devAtom.LIZStoreIdx, devOffsets,
                                                                kkrsz_ns, ispin, lsms.n_spin_pola, lsms.n_spin_cant,
                                                                iie, d.getBlkSizeTmatStore(), d.getTmatStoreLDim(),
                                                                (deviceDoubleComplex *)d.getDevTmatStore(), nrmat_ns,
                                                                (deviceDoubleComplex *)devBgij, (deviceDoubleComplex *)devM);
}

void buildKKRMatrixHip(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom,
                        DeviceStorage &devStorage, DeviceAtom &devAtom, int ispin,
                        int iie, Complex energy, Complex prel, Complex *devM)
{
  // decide between identical lmax and different lmax:

  // printf("buildKKRMatrixHip not finished yet!\n");
  // exit(1);
  
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

    // buildKKRMatrixLMaxIdenticalHip(lsms, local, atom, devStorage, devAtom, ispin,
    //                                iie, energy, prel, devM);
    buildKKRMatrixLMaxDifferentHip(lsms, local, atom, devStorage, devAtom, ispin,
                                    iie, energy, prel, devM);
  } else {
    // printf("lmax not identical in buildKKRMatrix\n");
    buildKKRMatrixLMaxDifferentHip(lsms, local, atom, devStorage, devAtom, ispin,
                                    iie, energy, prel, devM);
  }
}

