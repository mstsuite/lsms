/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

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

#ifdef ARCH_CUDA
#include "inversionTest_cuda.hpp"
#endif

extern "C"
{
    void block_inv_(Complex *a, Complex *vecs, int *lda, int *na, int *mp, int *ipvt, int *blk_sz, int *nblk, Complex *delta,
                  int *iwork, double *rwork, Complex *work1, int *alg,
                  int *idcol, int *iprint);

}

// #include "makegij_new.cpp"

void usage(const char *name)
{
  printf("usage: %s <matrix type> [options]\n",name);
  printf("  matrix type: 1: 1-tG and G, t Hilbert matrices, options: <block size> <num blocks>\n");
  printf("               2: 1-tG and G = -1, t = 1\n");
}

template <typename T>
void zeroMatrix(Matrix<T> &m)
{
  for(int i=0; i<m.n_row(); i++)
    for(int j=0; j<m.n_col(); j++)
      m(i,j) = 0.0;
}

template <typename T>
void unitMatrix(Matrix<T> &m)
{
  for(int i=0; i<m.n_row(); i++)
  {
    for(int j=0;j<m.n_col(); j++)
      m(i,j) = 0.0;
    m(i,i) = 1.0;
  }
}

Real matrixDistance(Matrix<Complex> &a, Matrix<Complex> &b)
{
  Real d;

  for(int i=0; i<a.n_col(); i++)
    for(int j=0; j<a.n_row(); j++)
      d += ((a(i,j)-b(i,j)) * std::conj(a(i,j)-b(i,j))).real();
  
  return std::sqrt(d);
}

void writeMatrix(Matrix<Complex> &a)
{
  for(int i=0; i<a.n_col(); i++)
    for(int j=0; j<a.n_row(); j++)
      printf("%d %d %f %f\n",i,j,a(i,j).real(), a(i,j).imag());
}

template <typename T>
void makeHilbertMatrix(Matrix<T> &m)
{
  for(int i=0; i<m.n_row(); i++)
    for(int j=0; j<m.n_col(); j++)
      m(i,j) = 1.0/(Complex(i+j+1));
}

template <typename T>
void zeroDiagonalBlocks(Matrix<T> &m, int blockSize)
{
  int n=m.n_row();
  for(int iBlock=0; iBlock<n/blockSize; iBlock++)
    for(int jBlock=0; jBlock<n/blockSize; jBlock++)
    {
      int ii=iBlock*blockSize;
      int jj=jBlock*blockSize;
      for(int i=0; i<std::min(blockSize, n-ii); i++)
        for(int j=0; j<std::min(blockSize, n-jj); j++)
          m(ii+i, jj+j) = 0.0;
    }
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

// given the  m-Matrix [(blockSize*numBlocks) x (blockSize*numBlocks)] and tMatrix[0] [blockSize x blockSize]
// calculate the t00-Matrix ast the blockSize x blockSize diagonal block of  (1 - tG)^-1 t
void solveTau00Reference(Matrix<Complex> &tau00, Matrix<Complex> &m, std::vector<Matrix<Complex> > &tMatrices, int blockSize, int numBlocks)
{
  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK
  Matrix<Complex> tau(blockSize * numBlocks, blockSize);
  zeroMatrix(tau);
  // copy t[0] into the top part of tau
  for(int i=0; i<blockSize; i++)
    for(int j=0; j<blockSize; j++)
      tau(i,j) = tMatrices[0](i,j);

  int n = blockSize * numBlocks;
  int ipiv[n];
  int info;
  LAPACK::zgesv_(&n, &blockSize, &m(0,0), &n, ipiv, &tau(0,0), &n, &info);

  // copy result into tau00
  for(int i=0; i<blockSize; i++)
    for(int j=0; j<blockSize; j++)
      tau00(i,j) = tau(i,j);
}

void solveTau00zgetrf(Matrix<Complex> &tau00, Matrix<Complex> &m, std::vector<Matrix<Complex> > &tMatrices, int blockSize, int numBlocks)
{
  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK
  Matrix<Complex> tau(blockSize * numBlocks, blockSize);
  
  zeroMatrix(tau);
  // copy t[0] into the top part of t
  for(int i=0; i<blockSize; i++)
    for(int j=0; j<blockSize; j++)
      tau(i,j) = tMatrices[0](i,j);

  int n = blockSize * numBlocks;
  int ipiv[n];
  Matrix<Complex> work(n,blockSize);
  std::vector<std::complex<float> > swork(n*(n+blockSize));
  std::vector<double> rwork(n);
  int info, iter;

  LAPACK::zgetrf_(&n, &n, &m(0,0), &n, &ipiv[0], &info);
  LAPACK::zgetrs_("N", &n, &blockSize, &m(0,0), &n, &ipiv[0], &tau(0,0), &n, &info);

  // copy result into tau00
  for(int i=0; i<blockSize; i++)
    for(int j=0; j<blockSize; j++)
      tau00(i,j) = tau(i,j);
}

#ifndef ARCH_IBM
void solveTau00zcgesv(Matrix<Complex> &tau00, Matrix<Complex> &m, std::vector<Matrix<Complex> > &tMatrices, int blockSize, int numBlocks)
{
  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK
  Matrix<Complex> tau(blockSize * numBlocks, blockSize);
  Matrix<Complex> t(blockSize * numBlocks, blockSize);
  
  zeroMatrix(tau);
  zeroMatrix(t);
  // copy t[0] into the top part of t
  for(int i=0; i<blockSize; i++)
    for(int j=0; j<blockSize; j++)
      t(i,j) = tMatrices[0](i,j);

  int n = blockSize * numBlocks;
  int ipiv[n];
  Matrix<Complex> work(n,blockSize);
  std::vector<std::complex<float> > swork(n*(n+blockSize));
  std::vector<double> rwork(n);
  int info, iter;
  LAPACK::zcgesv_(&n, &blockSize, &m(0,0), &n, ipiv, &t(0,0), &n, &tau(0,0), &n,
                  &work(0,0), &swork[0], &rwork[0],
                 &iter, &info);

  // copy result into tau00
  for(int i=0; i<blockSize; i++)
    for(int j=0; j<blockSize; j++)
      tau00(i,j) = tau(i,j);
}
#endif

void solveTau00zblocklu(Matrix<Complex> &tau00, Matrix<Complex> &m, std::vector<Matrix<Complex> > &tMatrices, int blockSize, int numBlocks)
{
  int nrmat_ns = blockSize * numBlocks;
  int ipvt[nrmat_ns];
  int info;
  int alg = 2, iprint = 0;
  Complex vecs[42]; // dummy, not used for alg=2
  int nblk = 3;
  Matrix<Complex> delta(blockSize, blockSize);
  int iwork[nrmat_ns];
  Real rwork[nrmat_ns];
  Complex work1[nrmat_ns];

  int blk_sz[1000];

  blk_sz[0]=blockSize;
  if(nblk==1)
    blk_sz[0]=nrmat_ns;
  else if(nblk==2)
    blk_sz[1]=nrmat_ns-blk_sz[0];
  else if(nblk>2)
  {
    int min_sz=(nrmat_ns-blk_sz[0])/(nblk-1);
    int rem=(nrmat_ns-blk_sz[0])%(nblk-1);
    int i=1;
    for(;i<=rem;i++)
      blk_sz[i]=min_sz+1;
    for(;i<nblk;i++)
      blk_sz[i]=min_sz;
  }

  int idcol[blk_sz[0]]; idcol[0]=0;

  // with m = [[A B][C D]], A: blk_sz[0] x blk_sz[0]
  // calculate the Schur complement m/D of m with A set to zero,
  // i.e. delta = B D^-1 C
  block_inv_(&m(0,0), vecs, &nrmat_ns, &nrmat_ns, &nrmat_ns, ipvt,
             blk_sz, &nblk, &delta(0,0),
             iwork, rwork, work1, &alg, idcol, &iprint);


  Matrix<Complex> wbig(blockSize, blockSize);
// setup unit matrix...............................................
// n.b. this is the top diagonal block of the kkr matrix m
//      i.e. 1 - t_0 G_00, with G_ii == 0 this is just the unit matrix

  unitMatrix(wbig);

// c     get 1-delta and put it in wbig

  for(int i=0; i<blockSize; i++)
    for(int j=0; j<blockSize; j++)
      wbig(i,j) -= delta(i,j);
//  c     ================================================================
// c     create tau00 => {[1-t*G]**(-1)}*t : for central site only.......
// c     ----------------------------------------------------------------

  LAPACK::zgetrf_(&blockSize, &blockSize, &wbig(0,0), &blockSize, ipvt, &info);

  for(int i=0; i<blockSize; i++)
    for(int j=0; j<blockSize; j++)
      tau00(i,j) = tMatrices[0](i,j);

  LAPACK::zgetrs_("N", &blockSize, &blockSize, &wbig(0,0), &blockSize, ipvt, &tau00(0,0), &blockSize, &info);

}

void block_inverse(Matrix<Complex> &a, int *blk_sz, int nblk, Matrix<Complex> &delta, int *ipvt, int *idcol);

void solveTau00zblocklu_cpp(Matrix<Complex> &tau00, Matrix<Complex> &m, std::vector<Matrix<Complex> > &tMatrices, int blockSize, int numBlocks)
{
  int nrmat_ns = blockSize * numBlocks;
  int ipvt[nrmat_ns];
  int info;
  int nblk = 3;
  Matrix<Complex> delta(blockSize, blockSize);

  int blk_sz[1000];

  blk_sz[0]=blockSize;
  if(nblk==1)
    blk_sz[0]=nrmat_ns;
  else if(nblk==2)
    blk_sz[1]=nrmat_ns-blk_sz[0];
  else if(nblk>2)
  {
    int min_sz=(nrmat_ns-blk_sz[0])/(nblk-1);
    int rem=(nrmat_ns-blk_sz[0])%(nblk-1);
    int i=1;
    for(;i<=rem;i++)
      blk_sz[i]=min_sz+1;
    for(;i<nblk;i++)
      blk_sz[i]=min_sz;
  }

  int idcol[blk_sz[0]]; idcol[0]=0;

  // with m = [[A B][C D]], A: blk_sz[0] x blk_sz[0]
  // calculate the Schur complement m/D of m with A set to zero,
  // i.e. delta = B D^-1 C
  block_inverse(m, blk_sz, nblk, delta, ipvt, idcol);

  Matrix<Complex> wbig(blockSize, blockSize);
// setup unit matrix...............................................
// n.b. this is the top diagonal block of the kkr matrix m
//      i.e. 1 - t_0 G_00, with G_ii == 0 this is just the unit matrix

  unitMatrix(wbig);

// c     get 1-delta and put it in wbig

  for(int i=0; i<blockSize; i++)
    for(int j=0; j<blockSize; j++)
      wbig(i,j) -= delta(i,j);
//  c     ================================================================
// c     create tau00 => {[1-t*G]**(-1)}*t : for central site only.......
// c     ----------------------------------------------------------------

  LAPACK::zgetrf_(&blockSize, &blockSize, &wbig(0,0), &blockSize, ipvt, &info);

  for(int i=0; i<blockSize; i++)
    for(int j=0; j<blockSize; j++)
      tau00(i,j) = tMatrices[0](i,j);

  LAPACK::zgetrs_("N", &blockSize, &blockSize, &wbig(0,0), &blockSize, ipvt, &tau00(0,0), &blockSize, &info);

}

int main(int argc, char *argv[])
{
  int matrixType, blockSize=18, numBlocks=113;
  bool printMatrices=false;

#ifdef ARCH_CUDA
  cublasHandle_t cublasHandle;
  initCuda(cublasHandle);
#endif
  
  printf("Test of inversion routine for LSMS\n");
  if(argc<2)
  {
    usage(argv[0]);
    return 0;
  }
  matrixType = atoi(argv[1]);
  if(argc>2)
  {
    blockSize = atoi(argv[2]);
  }
  if(argc>3)
  {
    numBlocks = atoi(argv[3]);
  }
  int n = blockSize*numBlocks;

  printf("Matrix type      : %6d\n",matrixType);
  printf("Block Size       : %6d\n",blockSize);
  printf("Number of Blocks : %6d\n",numBlocks);
  printf("Total Matrix Size: %6d\n",n);
  
  std::vector<Matrix<Complex> > tMatrices;
  Matrix<Complex> G0(n, n);
  
  tMatrices.resize(numBlocks);

#ifdef ARCH_CUDA
  DeviceData devData;
  allocDeviceData(devData, blockSize, numBlocks);
#endif

  if(matrixType == 1)
  {
    for(int i=0; i<numBlocks; i++)
    {
      tMatrices[i].resize(blockSize, blockSize);
      makeHilbertMatrix<Complex>(tMatrices[i]);
    }
    makeHilbertMatrix<Complex>(G0);
    zeroDiagonalBlocks(G0, blockSize);
  } else if(matrixType == 2) {
    for(int i=0; i<numBlocks; i++)
    {
      tMatrices[i].resize(blockSize, blockSize);
      unitMatrix<Complex>(tMatrices[i]);
    }
    zeroMatrix<Complex>(G0);
  } else {
    for(int i=0; i<numBlocks; i++)
    {
      tMatrices[i].resize(blockSize, blockSize);
      unitMatrix<Complex>(tMatrices[i]);
    }
    zeroMatrix(G0);
  }
  
  Matrix<Complex> m(n, n);
  
  makeType1Matrix(m, G0, tMatrices, blockSize, numBlocks);
  Matrix<Complex> tau00Reference(blockSize, blockSize);
  auto startTimeReference = std::chrono::system_clock::now();
  solveTau00Reference(tau00Reference, m, tMatrices, blockSize, numBlocks);
  auto endTimeReference = std::chrono::system_clock::now();
  std::chrono::duration<double> timeReference = endTimeReference - startTimeReference;

  makeType1Matrix(m, G0, tMatrices, blockSize, numBlocks);
  Matrix<Complex> tau00zgetrf(blockSize, blockSize);
  auto startTimeZgetrf = std::chrono::system_clock::now();
  solveTau00zgetrf(tau00zgetrf, m, tMatrices, blockSize, numBlocks);
  auto endTimeZgetrf = std::chrono::system_clock::now();
  std::chrono::duration<double> timeZgetrf = endTimeZgetrf - startTimeZgetrf;
  
  makeType1Matrix(m, G0, tMatrices, blockSize, numBlocks);
  Matrix<Complex> tau00zblocklu(blockSize, blockSize);
  auto startTimeZblocklu = std::chrono::system_clock::now();
  solveTau00zblocklu(tau00zblocklu, m, tMatrices, blockSize, numBlocks);
  auto endTimeZblocklu = std::chrono::system_clock::now();
  std::chrono::duration<double> timeZblocklu = endTimeZblocklu - startTimeZblocklu;

  makeType1Matrix(m, G0, tMatrices, blockSize, numBlocks);
  Matrix<Complex> tau00zblocklu_cpp(blockSize, blockSize);
  auto startTimeZblocklu_cpp = std::chrono::system_clock::now();
  solveTau00zblocklu_cpp(tau00zblocklu_cpp, m, tMatrices, blockSize, numBlocks);
  auto endTimeZblocklu_cpp = std::chrono::system_clock::now();
  std::chrono::duration<double> timeZblocklu_cpp = endTimeZblocklu_cpp - startTimeZblocklu_cpp;

#ifndef ARCH_IBM
  makeType1Matrix(m, G0, tMatrices, blockSize, numBlocks);
  Matrix<Complex> tau00zcgesv(blockSize, blockSize);
  auto startTimeZcgesv = std::chrono::system_clock::now();
  solveTau00zcgesv(tau00zcgesv, m, tMatrices, blockSize, numBlocks);
  auto endTimeZcgesv = std::chrono::system_clock::now();
  std::chrono::duration<double> timeZcgesv = endTimeZcgesv - startTimeZcgesv;
#endif

  Real d = matrixDistance(tau00Reference, tau00zblocklu);
  printf("d2 (t00Reference, tau00zblocklu) = %g\n", d);
  d = matrixDistance(tau00Reference, tau00zblocklu_cpp);
  printf("d2 (t00Reference, tau00zblocklu_cpp) = %g\n", d);
#ifndef ARCH_IBM
  d = matrixDistance(tau00Reference, tau00zcgesv);
  printf("d2 (t00Reference, tau00zcgesv) = %g\n", d);
#endif
  d = matrixDistance(tau00Reference, tau00zgetrf);
  printf("d2 (t00Reference, tau00zgetrf) = %g\n", d);
  
  printf("t(Reference) = %fsec\n",timeReference.count());
  printf("t(zblocklu)  = %fsec\n",timeZblocklu.count());
  printf("t(zblocklu_cpp)  = %fsec\n",timeZblocklu_cpp.count());
#ifndef ARCH_IBM
  printf("t(zcgesv)  = %fsec\n",timeZcgesv.count());
#endif
  printf("t(zgetrf)  = %fsec\n",timeZgetrf.count());

#ifdef ARCH_CUDA
  makeType1Matrix(m, G0, tMatrices, blockSize, numBlocks);
  Matrix<Complex> tau00zgetrf_cublas(blockSize, blockSize);
  auto startTimeZgetrf_cublas_transfer = std::chrono::system_clock::now();
  transferMatrixToGPU(devData.m, m);
  transferMatrixToGPU(devData.tMatrices[0], tMatrices[0]);
  auto startTimeZgetrf_cublas = std::chrono::system_clock::now();
  solveTau00zgetrf_cublas(cublasHandle, devData, tau00zgetrf_cublas, blockSize, numBlocks);
  auto endTimeZgetrf_cublas = std::chrono::system_clock::now();
  std::chrono::duration<double> timeZgetrf_cublas = endTimeZgetrf_cublas - startTimeZgetrf_cublas;
  std::chrono::duration<double> timeZgetrf_cublas_transfer = endTimeZgetrf_cublas - startTimeZgetrf_cublas_transfer;

  makeType1Matrix(m, G0, tMatrices, blockSize, numBlocks);
  Matrix<Complex> tau00zblocklu_cublas(blockSize, blockSize);
  auto startTimeZblocklu_cublas_transfer = std::chrono::system_clock::now();
  transferMatrixToGPU(devData.m, m);
  transferMatrixToGPU(devData.tMatrices[0], tMatrices[0]);
  auto startTimeZblocklu_cublas = std::chrono::system_clock::now();
  solveTau00zblocklu_cublas(cublasHandle, devData, tau00zblocklu_cublas, m, tMatrices, blockSize, numBlocks);
  auto endTimeZblocklu_cublas = std::chrono::system_clock::now();
  std::chrono::duration<double> timeZblocklu_cublas = endTimeZblocklu_cublas - startTimeZblocklu_cublas;
  std::chrono::duration<double> timeZblocklu_cublas_transfer = endTimeZblocklu_cublas - startTimeZblocklu_cublas_transfer;
  
  printf("\nCUDA and cuBLAS:\n");
  d = matrixDistance(tau00Reference, tau00zgetrf_cublas);
  printf("d2 (t00Reference, tau00zgetrf_cublas) = %g\n", d);
  printf("t(zgetrf_cublas)  = %fsec [%fsec]\n",timeZgetrf_cublas.count(), timeZgetrf_cublas_transfer.count());

  d = matrixDistance(tau00Reference, tau00zblocklu_cublas);
  printf("d2 (t00Reference, tau00zblocklu_cublas) = %g\n", d);
  printf("t(zblocklu_cublas)  = %fsec [%fsec]\n",timeZblocklu_cublas.count(), timeZblocklu_cublas_transfer.count());
#endif

  if(printMatrices)
  {
    printf("\ntau00Reference:\n"); writeMatrix(tau00Reference);
    printf("\ntau00zblocklu:\n"); writeMatrix(tau00zblocklu);
  }

#ifdef ARCH_CUDA
  freeDeviceData(devData);
  finalizeCuda(cublasHandle);
#endif
  
  return 0;
}
