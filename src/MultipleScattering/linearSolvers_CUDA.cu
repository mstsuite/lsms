/* -*- mode: C++; c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#include "linearSolvers.hpp"

#include <stdio.h>

#include "Complex.hpp"
#include "Matrix.hpp"
#include <vector>

#define IDX(i, j, lDim) (((j)*(lDim))+(i))

template <typename T>
void zeroMatrixCuda(T *devM, int lDim, int nCol)
{
//  for(int i=0; i<m.n_row(); i++)
//    for(int j=0; j<m.n_col(); j++)
//      m(i,j) = 0.0;
  cudaMemset(devM, 0, lDim*nCol*sizeof(T));
}

template <typename T>
__global__ void setDiagonalKernelCuda(T *devM, int lDim, int nCol, T val)
{
  int i=blockIdx.x;
  if(i<nCol)
  {
    devM[IDX(i, i, lDim)] = val;
  }
}

template <typename T>
__global__ void addDiagonalKernelCuda(T *devM, int lDim, int nCol, T val)
{
  int i=blockIdx.x;
  if(i<nCol)
  {
    devM[IDX(i, i, lDim)] = cuCadd(devM[IDX(i, i, lDim)], val);
  }
}

template <typename T>
void unitMatrixCuda(T *devM, int lDim, int nCol)
{
  zeroMatrixCuda(devM, lDim, nCol);
  setDiagonalKernelCuda<<<nCol,1>>>(devM, lDim, nCol, 1.0);
}

template <typename T>
__global__ void zeroDiagonalBlocksKernelCuda(T *devM, int lDim, int nCol, int blockSize)
{
  int iBlock = blockIdx.x;
  int jBlock = blockIdx.y;
  if(iBlock<nCol/blockSize)
    if(jBlock<nCol/blockSize)
    {
      int ii=iBlock*blockSize;
      int jj=jBlock*blockSize;
      for(int i=0; i<std::min(blockSize, nCol-ii); i++)
        for(int j=0; j<std::min(blockSize, nCol-jj); j++)
          devM[IDX(ii+i, jj+j, lDim)] = 0.0;
    }
}


void transferMatrixToGPUCuda(cuDoubleComplex *devM, Matrix<Complex> &m)
{
  cudaMemcpy(devM, &m(0,0), m.l_dim()*m.n_col()*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
}

void transferMatrixFromGPUCuda(Matrix<Complex> &m, cuDoubleComplex *devM)
{
  cudaMemcpy(&m(0,0), devM,  m.l_dim()*m.n_col()*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
}

__global__ void copyTMatrixToTauCuda(cuDoubleComplex *tau, cuDoubleComplex *t, int kkrsz, int nrmat)
{
  int i = blockIdx.x;
  if(i < kkrsz)
  {
    for(int j=0; j<kkrsz; j++)
      tau[IDX(i,j,nrmat)] = t[IDX(i,j,kkrsz)];
  }
}

__global__ void copyTauToTau00Cuda(cuDoubleComplex *tau00, cuDoubleComplex *tau, int kkrsz, int nrmat)
{
  int i = blockIdx.x;
  if(i < kkrsz)
  {
    for(int j=0; j<kkrsz; j++)
      tau00[IDX(i,j,kkrsz)] = tau[IDX(i,j,nrmat)];
  }
}

void solveTau00zgetrf_cublas(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom,
                             Complex &tMatrix, Complex *devM,
                             Matrix<Complex> &tau00)
{
  cublasHandle_t cublasHandle = DeviceStorage::getCublasHandle()
  int nrmat_ns = lsms.n_spin_cant*atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block
  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK
  cuDoubleComplex *Aarray[1], *Barray[1];

  cuDoubleComplex *devTau = d.getDevTau();
  cuDoubleComplex *devTau00 = d.getDevTau00();
  // printf("zero Matrix\n");
  zeroMatrixCuda(devTau, nrmat_ns, kkrzs_ns);
  // printf("copyTMatrixToTau\n");
  copyTMatrixToTauCuda<<<blockSize,1>>>(devTau, tMatrix, kkrsz_ns, nrmat_ns);

  Barray[0] = devTau;

  Aarray[0] = devM;
  
  int *ipivArray=d.getDevIpiv;
  int *infoArray=d.info;
  int info;

  // printf("cublasZgetrfBatched\n");
  cublasZgetrfBatched(cublasHandle, nrmat_ns, Aarray, nrmat_ns, ipivArray, infoArray, 1);
  // printf("cublasZgetrsBatched\n");

  cublasZgetrsBatched(cublasHandle, CUBLAS_OP_N, nrmat_ns, kkrzs_ns, Aarray, nrmat_ns, ipivArray,
                      Barray, nrmat_ns, &info, 1);

  // copy result into tau00
  // printf("copyTauToTau00\n");
  copyTauToTau00Cuda<<<blockSize,1>>>(devTau00, devTau, kkrzs_ns, nrmat_ns);
  // printf("transferMatrixFromGPU\n");
  transferMatrixFromGPUCuda(tau00, devTau00);
}

#ifndef ARCH_IBM
void solveTau00zzgesv_cusolver(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom,
                               Complex *tMatrix, Complex *devM, Matrix<Complex> &tau00)
{
  cusolverDnHandle_t cusolverDnHandle = DeviceStorage::getCusolverDnHandle();
  int nrmat_ns = lsms.n_spin_cant*atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block
  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK

  cuDoubleComplex *devTau = d.getDevTau();
  cuDoubleComplex *devTau00 = d.getDevTau00();
  
  zeroMatrixCuda(devTau, blockSize*numBlocks, kkrsz_ns);
  zeroMatrixCuda(deviceData.t, blockSize*numBlocks, kkrsz_ns);
  copyTMatrixToTauCuda<<<blockSize,1>>>(deviceData.t, tMatrix], kkrsz_ns, nrmat_ns);

  int info, iter;

  cusolverStatus_t status = cusolverDnZZgesv(cusolverDnHandle, nrmat_ns, kkrsz_ns,
                                             devM, nrmat_ns, deviceData.ipiv, deviceData.t, nrmat_ns, deviceData.tau, nrmat_ns,
    deviceData.work, deviceData.workBytes, &iter, deviceData.info);

  if(status!=CUSOLVER_STATUS_SUCCESS)
  {
    printf("cusolverDnZZgesv returned %d\n",status);
  }

  copyTauToTau00Cuda<<<blockSize,1>>>(deviceData.tau00, deviceData.tau, kkrsz_ns, nrmat_ns);
  transferMatrixFromGPU(tau00, deviceData.tau00);
}
#endif

void solveTau00zgetrf_cusolver(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom,
                               Complex *tMatrix, Complex *devM, Matrix<Complex> &tau00)
{
  cusolverDnHandle_t cusolverDnHandle = DeviceStorage::getCusolverDnHandle();
  int nrmat_ns = lsms.n_spin_cant*atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block
  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK

  zeroMatrixCuda(deviceData.tau, nrmat_ns, kkrsz_ns);
  copyTMatrixToTau<<<blockSize,1>>>(deviceData.tau, tMatrix, kkrsz_ns, nrmat_ns);

  cusolverDnZgetrf(cusolverDnHandle, nrmat_ns, nrmat_ns, 
                   devM, nrmat_ns,
           (cuDoubleComplex *)deviceData.work,
           deviceData.ipiv,
           deviceData.info );

  cusolverDnZgetrs(cusolverDnHandle, CUBLAS_OP_N, nrmat_ns, kkrsz_ns,
      devM, nrmat_ns, deviceData.ipiv, deviceData.tau, nrmat_ns, deviceData.info);

  // copy result into tau00
  copyTauToTau00<<<blockSize,1>>>(deviceData.tau00, deviceData.tau, kkrsz_ns, nrmat_ns);
  transferMatrixFromGPU(tau00, deviceData.tau00);
}
