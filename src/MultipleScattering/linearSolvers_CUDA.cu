/* -*- mode: C++; c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#include "linearSolvers.hpp"

#include <stdio.h>

#include "Complex.hpp"
#include "Matrix.hpp"
#include <vector>

#include "Accelerator/DeviceStorage.hpp"
#include "Accelerator/deviceCheckError.hpp"
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cublas_v2.h>
#include <cusolverDn.h>

/*
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
  int i=blockIdx.x*blockDim.x + threadIdx.x;
  if(i<nCol)
  {
    devM[IDX(i, i, lDim)] = val;
  }
}

template <typename T>
__global__ void addDiagonalKernelCuda(T *devM, int lDim, int nCol, T val)
{
  int i=blockIdx.x*blockDim.x + threadIdx.x;
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

*/



template <typename T>
__global__ void zeroDiagonalBlocksKernelCuda(T *devM, int lDim, int nCol, int blockSize)
{
  int iBlock = blockIdx.x*blockDim.x + threadIdx.x;
  int jBlock = blockIdx.y*blockDim.y + threadIdx.y;
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

void transferT0MatrixToGPUCuda(Complex *devT0, LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie)
{
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz;
  cudaMemcpy(devT0, &local.tmatStore(iie*local.blkSizeTmatStore,atom.LIZStoreIdx[0]),
             kkrsz_ns*kkrsz_ns*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
}

void transferMatrixToGPUCuda(Complex *devM, Matrix<Complex> &m)
{
  cudaMemcpy(devM, &m(0,0), m.l_dim()*m.n_col()*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
}

void transferMatrixFromGPUCuda(Matrix<Complex> &m, cuDoubleComplex *devM)
{
  cudaMemcpy(&m(0,0), devM,  m.l_dim()*m.n_col()*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
}

__global__ void copyTMatrixToTauCuda(cuDoubleComplex *tau, cuDoubleComplex *t, int kkrsz, int nrmat)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i < kkrsz)
  {
    for(int j=0; j<kkrsz; j++)
      tau[IDX(i,j,nrmat)] = t[IDX(i,j,kkrsz)];
  }
}

__global__ void copyTauToTau00Cuda(cuDoubleComplex *tau00, cuDoubleComplex *tau, int kkrsz, int nrmat)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i < kkrsz)
  {
    for(int j=0; j<kkrsz; j++)
      tau00[IDX(i,j,kkrsz)] = tau[IDX(i,j,nrmat)];
  }
}

void solveTau00zgetrf_cublas(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom,
                             Complex *tMatrix, Complex *devM,
                             Matrix<Complex> &tau00)
{
  cublasHandle_t cublasHandle = DeviceStorage::getCublasHandle();
  int nrmat_ns = lsms.n_spin_cant*atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block
  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK
  cuDoubleComplex *Aarray[1], *Barray[1];

  cuDoubleComplex *devTau = (cuDoubleComplex *)d.getDevTau();
  cuDoubleComplex *devTau00 = (cuDoubleComplex *)d.getDevTau00();
  // printf("zero Matrix\n");
  zeroMatrixCuda(devTau, nrmat_ns, kkrsz_ns);
  deviceCheckError();
  // printf("copyTMatrixToTau\n");
  copyTMatrixToTauCuda<<<kkrsz_ns,1>>>(devTau, (cuDoubleComplex *)tMatrix, kkrsz_ns, nrmat_ns);
  deviceCheckError();
  
  Barray[0] = devTau;

  Aarray[0] = (cuDoubleComplex *)devM;
  
  int *ipivArray=d.getDevIpvt();
  int *infoArray = d.getDevInfo();
  int info;

  // printf("cublasZgetrfBatched\n");
  cublasCheckError(cublasZgetrfBatched(cublasHandle, nrmat_ns, Aarray, nrmat_ns, ipivArray, infoArray, 1));
  // printf("cublasZgetrsBatched\n");

  cublasCheckError(cublasZgetrsBatched(cublasHandle, CUBLAS_OP_N, nrmat_ns, kkrsz_ns, Aarray, nrmat_ns, ipivArray,
                                 Barray, nrmat_ns, &info, 1));

  // copy result into tau00
  // printf("copyTauToTau00\n");
  copyTauToTau00Cuda<<<kkrsz_ns,1>>>(devTau00, devTau, kkrsz_ns, nrmat_ns);
  deviceCheckError();
  // printf("transferMatrixFromGPU\n");
  transferMatrixFromGPUCuda(tau00, devTau00);
  deviceCheckError();
}

#ifndef ARCH_IBM
void solveTau00zzgesv_cusolver(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom,
                               Complex *tMatrix, Complex *devM, Matrix<Complex> &tau00)
{
  cusolverDnHandle_t cusolverDnHandle = DeviceStorage::getCusolverDnHandle();
  int nrmat_ns = lsms.n_spin_cant*atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block
  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK

  cuDoubleComplex *devTau = (cuDoubleComplex *)d.getDevTau();
  cuDoubleComplex *devTau00 = (cuDoubleComplex *)d.getDevTau00();
  cuDoubleComplex *devWork = (cuDoubleComplex *)d.getDevWork();

  cuDoubleComplex *devT = (cuDoubleComplex *)d.getDevT();

  int *devIpiv = d.getDevIpvt();
  int devInfo[1]; // d.getDevInfo();
  
  zeroMatrixCuda(devTau, nrmat_ns, kkrsz_ns);
  zeroMatrixCuda(devT, nrmat_ns, kkrsz_ns);
  copyTMatrixToTauCuda<<<kkrsz_ns,1>>>(devT, (cuDoubleComplex *)tMatrix, kkrsz_ns, nrmat_ns);

  int iter;

  cusolverStatus_t status = cusolverDnZZgesv(cusolverDnHandle, nrmat_ns, kkrsz_ns,
                                             (cuDoubleComplex *)devM, nrmat_ns, devIpiv, devT, nrmat_ns, devTau, nrmat_ns,
                                             devWork, d.getDevWorkBytes(), &iter, devInfo);

  if(status!=CUSOLVER_STATUS_SUCCESS)
  {
    printf("cusolverDnZZgesv returned %d\n",status);
  }

  copyTauToTau00Cuda<<<kkrsz_ns,1>>>(devTau00, devTau, kkrsz_ns, nrmat_ns);
  transferMatrixFromGPUCuda(tau00, devTau00);
}
#endif

void solveTau00zgetrf_cusolver(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom,
                               Complex *tMatrix, Complex *devM, Matrix<Complex> &tau00)
{
  cusolverDnHandle_t cusolverDnHandle = DeviceStorage::getCusolverDnHandle();
  int nrmat_ns = lsms.n_spin_cant*atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block
  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK
  cuDoubleComplex *devTau = (cuDoubleComplex *)d.getDevTau();
  cuDoubleComplex *devTau00 = (cuDoubleComplex *)d.getDevTau00();
  cuDoubleComplex *devWork = (cuDoubleComplex *)d.getDevWork();

  int *devIpiv=d.getDevIpvt();
  int *devInfo = d.getDevInfo();

  zeroMatrixCuda(devTau, nrmat_ns, kkrsz_ns);
  deviceCheckError();
  copyTMatrixToTauCuda<<<kkrsz_ns,1>>>(devTau, (cuDoubleComplex *)tMatrix, kkrsz_ns, nrmat_ns);
  deviceCheckError();
  
  cusolverCheckError(cusolverDnZgetrf(cusolverDnHandle, nrmat_ns, nrmat_ns, 
                                      (cuDoubleComplex *)devM, nrmat_ns, devWork, devIpiv,
                                      devInfo ));

  cusolverCheckError(cusolverDnZgetrs(cusolverDnHandle, CUBLAS_OP_N, nrmat_ns, kkrsz_ns,
                                      (cuDoubleComplex *)devM, nrmat_ns, devIpiv, devTau, nrmat_ns, devInfo));
  
  // copy result into tau00
  copyTauToTau00Cuda<<<kkrsz_ns,1>>>(devTau00, devTau, kkrsz_ns, nrmat_ns);
  deviceCheckError();
  transferMatrixFromGPUCuda(tau00, devTau00);
  deviceCheckError();
}

#ifdef USE_XGETRF
void solveTau00Xgetrf_cusolver(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom,
                               Complex *tMatrix, Complex *devM, Matrix<Complex> &tau00)
{
  cusolverDnHandle_t cusolverDnHandle = DeviceStorage::getCusolverDnHandle();
  cusolverDnParams_t cusolverDnParams = DeviceStorage::getCusolverParams();
  int nrmat_ns = lsms.n_spin_cant*atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block
  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK
  cuDoubleComplex *devTau = (cuDoubleComplex *)d.getDevTau();
  cuDoubleComplex *devTau00 = (cuDoubleComplex *)d.getDevTau00();
  void *devWork = d.getDevWork();
  size_t devWorkBytes = d.getDevWorkBytes();
  int64_t *devIpiv=d.getDevIpvt64();
  void *hostWork = d.getHostWork();
  size_t hostWorkBytes = d.getHostWorkBytes();
  int *devInfo = d.getDevInfo();

  zeroMatrixCuda(devTau, nrmat_ns, kkrsz_ns);
  deviceCheckError();
  copyTMatrixToTauCuda<<<kkrsz_ns,1>>>(devTau, (cuDoubleComplex *)tMatrix, kkrsz_ns, nrmat_ns);
  deviceCheckError();

  cusolverCheckError(cusolverDnXgetrf(cusolverDnHandle,
                                      cusolverDnParams,
                                      (int64_t)nrmat_ns,
                                      (int64_t)nrmat_ns,
                                      CUDA_C_64F,
                                      (cuDoubleComplex *)devM,
                                      (int64_t)nrmat_ns,
                                      devIpiv,
                                      CUDA_C_64F,
                                      devWork,
                                      devWorkBytes,
                                      hostWork,
                                      hostWorkBytes,
                                      devInfo));

  cusolverCheckError(cusolverDnXgetrs(cusolverDnHandle,
                                      cusolverDnParams,
                                      CUBLAS_OP_N,
                                      (int64_t)nrmat_ns,
                                      (int64_t)kkrsz_ns,
                                      CUDA_C_64F,
                                      (cuDoubleComplex *)devM,
                                      (int64_t)nrmat_ns,
                                      devIpiv,
                                      CUDA_C_64F,
                                      devTau,
                                      (int64_t)nrmat_ns,
                                      devInfo));
  
  // copy result into tau00
  copyTauToTau00Cuda<<<kkrsz_ns,1>>>(devTau00, devTau, kkrsz_ns, nrmat_ns);
  deviceCheckError();
  transferMatrixFromGPUCuda(tau00, devTau00);
  deviceCheckError();
}
#endif

#ifdef USE_IRSXGESV
void solveTau00IRSXgesv_cusolver(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom,
                                 Complex *tMatrix, Complex *devM, Matrix<Complex> &tau00)
{
  cusolverDnHandle_t cusolverDnHandle = DeviceStorage::getCusolverDnHandle();
  cusolverDnIRSParams_t cusolverDnIRSParams = DeviceStorage::getCusolverIRSParams();
  cusolverDnIRSInfos_t cusolverDnIRSInfo = DeviceStorage::getCusolverIRSInfo();
  int nrmat_ns = lsms.n_spin_cant*atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block
  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK
  cuDoubleComplex *devTau = (cuDoubleComplex *)d.getDevTau();
  cuDoubleComplex *devTau00 = (cuDoubleComplex *)d.getDevTau00();
  cuDoubleComplex *devWork = (cuDoubleComplex *)d.getDevWork();
  cuDoubleComplex *devX = (cuDoubleComplex *)d.getDevX();
  size_t devWorkBytes = d.getDevWorkBytes();
  int *devInfo = d.getDevInfo();

  zeroMatrixCuda(devTau, nrmat_ns, kkrsz_ns);
  deviceCheckError();
  copyTMatrixToTauCuda<<<kkrsz_ns,1>>>(devTau, (cuDoubleComplex *)tMatrix, kkrsz_ns, nrmat_ns);
  deviceCheckError();

  cusolverDnIRSInfos_t info;
  int niters;
  cusolverCheckError(cusolverDnIRSXgesv(cusolverDnHandle,
                                        cusolverDnIRSParams,
                                        cusolverDnIRSInfo,
                                        nrmat_ns,
                                        kkrsz_ns,
                                        (cuDoubleComplex *)devM,
                                        nrmat_ns,
                                        devTau,
                                        nrmat_ns,
                                        devX,
                                        nrmat_ns,
                                        devWork,
                                        devWorkBytes,
                                        &niters,
                                        devInfo));

  // copy result into tau00
  copyTauToTau00Cuda<<<kkrsz_ns,1>>>(devTau00, devX, kkrsz_ns, nrmat_ns);
  deviceCheckError();
  transferMatrixFromGPUCuda(tau00, devTau00);
  deviceCheckError();
}
#endif
