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

void allocDeviceData(DeviceData &d, int blockSize, int numBlocks)
{
  int n=blockSize*numBlocks;
  cudaMalloc((void**)&d.m, n*n*sizeof(cuDoubleComplex));
  cudaMalloc((void**)&d.G0, n*n*sizeof(cuDoubleComplex));
  cudaMalloc((void**)&d.tau, n*blockSize*sizeof(cuDoubleComplex));
  cudaMalloc((void**)&d.tau00, blockSize*blockSize*sizeof(cuDoubleComplex));
  cudaMalloc((void**)&d.tMatrixStore, blockSize*blockSize*numBlocks*sizeof(cuDoubleComplex));
  tMatrices.resize(numBlocks);
  for(int i=0; i<numBlocks; i++)
    tMatrices[i] = &tMatrixStore[blockSize * blockSize * i];
  cudaMalloc((void**)&ipiv, n*sizeof(int));
  cudaMalloc((void**)&info, sizeof(int));
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

void usage(const char *name)
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
  cudaMemset(devM, 0, lDim*nCol*sizeof<T>);
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
void unitMatrixCuda(T *devM, int lDim, int nCol)
{
  zeroMatrixCuda(T, lDim, nCol);
  setDiagonalKernel<<<nCol,1>>>(devM, lDim, nCol, 1.0);
}

Real matrixDistance(Matrix<Complex> &a, Matrix<Complex> &b)
{
  Real d;

  for(int i=0; i<a.n_col(); i++)
    for(int j=0; j<a.n_row(); j++)
      d += ((a(i,j)-b(i,j)) * std::conj(a(i,j)-b(i,j))).real();
  
  return std::sqrt(d);
}


template <typename T>
__global__ void makeHilbertMatrixKernel(T *devM, int lDim, int nCol)
{
  int i = blockIdx.x;
  if(i<lDim)
  {
    for(int j=0; j<nCol; j++)
      devM[IDX(i,j, lDim)] = 1.0/(Complex(i+j+1));
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
          m[IDX(ii+i, jj+j, lDim)] = 0.0;
    }
}

template <typename T>
void zeroDiagonalBlocksCuda(T *devM, int lDim, int nCol, int blockSize)
{
  zeroDiagonalBlocksKernel<<<nCol/blockSize,nCol/blockSize>>>(devM, lDim, nCol, blockSize);
}

// type 1 matrix:
//
void makeType1Matrix(Matrix<Complex> &m, Matrix<Complex> &G0, std::vector<Matrix<Complex> > &tMatrices, int blockSize, int numBlocks)
{
  int n = m.n_row();
  Complex mone = -1.0;
  Complex zero = 0.0;
  // unitMatrix(m);
  // loop over the blocks to build -tG
  // m_ij <- t_i G_ij
  for(int iBlock=0; iBlock<numBlocks; iBlock++)
    for(int jBlock=0; jBlock<numBlocks; jBlock++)
    {
      // ZGEMM(TRANSA,TRANSB,M,N,K,ALPHA, A,LDA, B,LDB, BETA,C,LDC)
      // C := alpha*op( A )*op( B ) + beta*C,
      BLAS::zgemm_("N","N",&blockSize,&blockSize,&blockSize,&mone,
             &tMatrices[iBlock](0,0), &blockSize,
             &G0(iBlock*blockSize,jBlock*blockSize), &n,
             &zero, &m(iBlock*blockSize,jBlock*blockSize), &n);
    }
  // add unit matrix
  for(int i=0; i<m.n_row(); i++)
    m(i,i) = 1.0 + m(i,i);
}

void makeType1MatrixGPU(cublasHandle_t cublasHandle, DeviceData &d, int blockSize, int numBlocks)
{
  int n=blockSize*numBlocks;
  cuDoubleComplex mone = -1.0;
  cuDoubleComplex zero = 0.0;
  std::vector<cuDoubleComplex *> ts(numBlocks*numBlocks);
  std::vector<cuDoubleComplex *> G0s(numBlocks*numBlocks);
  std::vector<cuDoubleComplex *> ms(numBlocks*numBlocks);

  for(int iBlock=0; iBlock<numBlocks; iBlock++)
    for(int jBlock=0; jBlock<numBlocks; jBlock++)
    {
      ts[iBlock + jBlock*numBlocks] = &d.tMatrices[iBlock];
      G0s[iBlock + jBlock*numBlocks] = &d.G0[IDX(iBlock*blockSize,jBlock*blockSize,n)];
      m[iBlock + jBlock*numBlocks] = &d.m[IDX(iBlock*blockSize,jBlock*blockSize,n)];
    }

  cublasZgemmBatched(cublasHandle, CUBLAS_OP_N, CUBLAS_OP_N, blockSize, blockSize, blockSize,
                                  &mone,
                                  const cuDoubleComplex *Aarray[], int lda,
                                  const cuDoubleComplex *Barray[], int ldb,
                                  const cuDoubleComplex *beta,
                                  cuDoubleComplex *Carray[], int ldc, 
                                  int batchCount)
}

void transferMatrixToGPU(Complex *devM, Matrix<Complex> &m)
{
  cudaMemcpy(devM, &m(0,0), m.l_dim()*n_col()*sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
}

void transferMatrixFromGPU(Matrix<Complex> &m, Complex *devM)
{
  cudaMemcpy(&m(0,0), devM,  m.l_dim()*n_col()*sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
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
  
  zeroMatrixGPU(d.tau);
  copyTMatrixToTau<<<blockSize,1>>>(d.tau, d.tMatrices[0], blockSize, numBlocks);

  Barray[0] = d.tau;
  
  int n = blockSize * numBlocks;
  int *ipivArray, *infoArray;
  cudaMalloc((void**)&ipiv, n * sizeof(int));
  cudaMalloc((void**)&infoArray, 1 * sizeof(int));

  cublasZgetrfBatched(cublasHandle, n, Aarray, n, ipivArray, infoArray, 1);  
  cublasZgetrsBatched(cublasHandle, CUBLAS_OP_N, n, blockSize, Aarray, n, ipivArray,
                      Barray, n, infoArray, 1);

  // copy result into tau00
  copyTauToTau00<<<blockSize,1>>>(d.tau00, d.tau, blockSize, numBlocks);
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


