/* -*- mode: C++; c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

// test the inverion algorithm for multiple scattering codes for the solution of
// tau = (1 - tG)^-1 t
// where t is a block diagonal matrix
// note that the diagonal blocks G_ii == 0

#include <stdio.h>

#include "Complex.hpp"
#include "Matrix.hpp"
#include <vector>
#include <chrono>
#include <ctime>

#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cublas_v2.h>

#include "inversionTest_cuda.hpp"

/*
class DeviceData
{
public:
  std::vector<cuDoubleComplex *> tMatrices;
  cuDoubleComplex *tMatrixStore;
  cuDoubleComplex *tau;
  cuDoubleComplex *tau00;
  cuDoubleComplex *m;
  cuDoubleComplex *G0;
  int *ipiv;
  int *info;
};
*/

void allocDeviceData(DeviceData &d, int blockSize, int numBlocks)
{
  int n=blockSize*numBlocks;
  cudaMalloc((void**)&d.m, n*n*sizeof(cuDoubleComplex));
  cudaMalloc((void**)&d.G0, n*n*sizeof(cuDoubleComplex));
  cudaMalloc((void**)&d.tau, n*blockSize*sizeof(cuDoubleComplex));
  cudaMalloc((void**)&d.tau00, blockSize*blockSize*sizeof(cuDoubleComplex));
  cudaMalloc((void**)&d.tMatrixStore, blockSize*blockSize*numBlocks*sizeof(cuDoubleComplex));
  d.tMatrices.resize(numBlocks);
  for(int i=0; i<numBlocks; i++)
    d.tMatrices[i] = &d.tMatrixStore[blockSize * blockSize * i];
  cudaMalloc((void**)&d.ipiv, n*sizeof(int));
  cudaMalloc((void**)&d.info, sizeof(int));
}

void freeDeviceData(DeviceData &d)
{
  cudaFree(d.tau);
  cudaFree(d.tau00);
  cudaFree(d.m);
  cudaFree(d.G0);
  cudaFree(d.tMatrixStore);
  cudaFree(d.ipiv);
  cudaFree(d.info);
}

// #include "makegij_new.cpp"

void usage_cuda(const char *name)
{
  printf("usage: %s <matrix type> [options]\n",name);
  printf("  matrix type: 1: 1-tG and G, t Hilbert matrices, options: <block size> <num blocks>\n");
  printf("               2: 1-tG and G = -1, t = 1\n");
}

#define IDX(i, j, lDim) (((j)*(lDim))+i)

template <typename T>
void zeroMatrixCuda(T *devM, int lDim, int nCol)
{
//  for(int i=0; i<m.n_row(); i++)
//    for(int j=0; j<m.n_col(); j++)
//      m(i,j) = 0.0;
  cudaMemset(devM, 0, lDim*nCol*sizeof(T));
}

template <typename T>
__global__ void setDiagonalKernel(T *devM, int lDim, int nCol, T val)
{
  int i=blockIdx.x;
  if(i<nCol)
  {
    devM[IDX(i, i, lDim)] = val;
  }
}

template <typename T>
__global__ void addDiagonalKernel(T *devM, int lDim, int nCol, T val)
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
  setDiagonalKernel<<<nCol,1>>>(devM, lDim, nCol, 1.0);
}

/*
Real matrixDistance(Matrix<Complex> &a, Matrix<Complex> &b)
{
  Real d;

  for(int i=0; i<a.n_col(); i++)
    for(int j=0; j<a.n_row(); j++)
      d += ((a(i,j)-b(i,j)) * std::conj(a(i,j)-b(i,j))).real();
  
  return std::sqrt(d);
}
*/

template <typename T>
__global__ void makeHilbertMatrixKernel(T *devM, int lDim, int nCol)
{
  int i = blockIdx.x;
  if(i<lDim)
  {
    for(int j=0; j<nCol; j++)
      devM[IDX(i,j, lDim)] = 1.0/(Complex(i+j+1));
  }
}

template <typename T>
__global__ void zeroDiagonalBlocksKernel(T *devM, int lDim, int nCol, int blockSize)
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

template <typename T>
void zeroDiagonalBlocksCuda(T *devM, int lDim, int nCol, int blockSize)
{
  zeroDiagonalBlocksKernel<<<nCol/blockSize,nCol/blockSize>>>(devM, lDim, nCol, blockSize);
}

void makeType1MatrixGPU(cublasHandle_t cublasHandle, DeviceData &d, int blockSize, int numBlocks)
{
  int n=blockSize*numBlocks;
  cuDoubleComplex one = {1.0, 0.0};
  cuDoubleComplex mone = {-1.0, 0.0};
  cuDoubleComplex zero = {0.0, 0.0};
  std::vector<cuDoubleComplex *> ts(numBlocks*numBlocks);
  std::vector<cuDoubleComplex *> G0s(numBlocks*numBlocks);
  std::vector<cuDoubleComplex *> ms(numBlocks*numBlocks);

  for(int iBlock=0; iBlock<numBlocks; iBlock++)
  {
    for(int jBlock=0; jBlock<numBlocks; jBlock++)
    {
      ts[iBlock + jBlock*numBlocks] = d.tMatrices[iBlock];
      G0s[iBlock + jBlock*numBlocks] = &d.G0[IDX(iBlock*blockSize,jBlock*blockSize,n)];
      ms[iBlock + jBlock*numBlocks] = &d.m[IDX(iBlock*blockSize,jBlock*blockSize,n)];
    }
  }

  cublasZgemmBatched(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, blockSize, blockSize, blockSize, &mone,
                                  &ts[0], blockSize,
                                  &G0s[0], n,
                                  &zero,
                                  &ms[0], n, 
                                  numBlocks*numBlocks);
  addDiagonalKernel<<<n,1>>>(d.m, n, n, one);
}

void transferMatrixToGPU(cuDoubleComplex *devM, Matrix<Complex> &m)
{
  cudaMemcpy(devM, &m(0,0), m.l_dim()*m.n_col()*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
}

void transferMatrixFromGPU(Matrix<Complex> &m, cuDoubleComplex *devM)
{
  cudaMemcpy(&m(0,0), devM,  m.l_dim()*m.n_col()*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
}

__global__ void copyTMatrixToTau(cuDoubleComplex *tau, cuDoubleComplex *t, int blockSize, int numBlocks)
{
  int i = blockIdx.x;
  int n = blockSize*numBlocks;
  if(i < blockSize)
  {
    for(int j=0; j<blockSize; j++)
      tau[IDX(i,j,n)] = t[IDX(i,j,blockSize)];
  }
}

__global__ void copyTauToTau00(cuDoubleComplex *tau00, cuDoubleComplex *tau, int blockSize, int numBlocks)
{
  int i = blockIdx.x;
  int n = blockSize*numBlocks;
  if(i < blockSize)
  {
    for(int j=0; j<blockSize; j++)
      tau00[IDX(i,j,blockSize)] = tau[IDX(i,j,n)];
  }
}

void solveTau00zgetrf_cublas(cublasHandle_t cublasHandle, DeviceData &d,
                             Matrix<Complex> &tau00, int blockSize, int numBlocks)
{
  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK
  cuDoubleComplex *Aarray[1], *Barray[1];
 
  // printf("zero Matrix\n"); 
  zeroMatrixCuda(d.tau, blockSize*numBlocks, blockSize);
  // printf("copyTMatrixToTau\n");
  copyTMatrixToTau<<<blockSize,1>>>(d.tau, d.tMatrices[0], blockSize, numBlocks);

  Barray[0] = d.tau;

  Aarray[0] = d.m;
  
  int n = blockSize * numBlocks;
  int *ipivArray=d.ipiv;
  int *infoArray=d.info;
  int info;

  // printf("cublasZgetrfBatched\n");
  cublasZgetrfBatched(cublasHandle, n, Aarray, n, ipivArray, infoArray, 1);
  // printf("cublasZgetrsBatched\n");

  cublasZgetrsBatched(cublasHandle, CUBLAS_OP_N, n, blockSize, Aarray, n, ipivArray,
                      Barray, n, &info, 1);

  // copy result into tau00
  // printf("copyTauToTau00\n");
  copyTauToTau00<<<blockSize,1>>>(d.tau00, d.tau, blockSize, numBlocks);
  // printf("transferMatrixFromGPU\n");
  transferMatrixFromGPU(tau00, d.tau00);
}

void initCuda(cublasHandle_t &cublasHandle)
{
  cublasCreate(&cublasHandle);
}

void finalizeCuda(cublasHandle_t &cublasHandle)
{
  cublasDestroy(cublasHandle);
}


