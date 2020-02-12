/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

// test the inverion algorithm for multiple scattering codes for the solution of
// tau = (1 - tG)^-1 t
// where t is a block diagonal matrix
// note that the diagonal blocks G_ii == 0

#include <stdio.h>

#include "Complex.hpp"
#include "Matrix.hpp"
#include <vector>

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
  
  block_inv_(&m(0,0), vecs, &nrmat_ns, &nrmat_ns, &nrmat_ns, ipvt,
             blk_sz, &nblk, &delta(0,0),
             iwork, rwork, work1, &alg, idcol, &iprint);


  Matrix<Complex> wbig(blockSize, blockSize);
//   c     setup unit matrix...............................................
// c     ----------------------------------------------------------------
//       call cmtruni(wbig,kkrsz_ns)
  unitMatrix(wbig);
//  zeroMatrix(wbig);
// c     ----------------------------------------------------------------
// c     get 1-delta and put it in wbig
//      call zaxpy(mtxsize,cmone,delta,1,wbig,1)
  for(int i=0; i<blockSize; i++)
    for(int j=0; j<blockSize; j++)
      wbig(i,j) -= delta(i,j);
//  c     ================================================================
// c     create tau00 => {[1-t*G]**(-1)}*t : for central site only.......
// c     ----------------------------------------------------------------
//      call zgetrf(kkrsz_ns,kkrsz_ns,wbig,kkrsz_ns,ipvt,info)
  LAPACK::zgetrf_(&blockSize, &blockSize, &wbig(0,0), &blockSize, ipvt, &info);
//      call zcopy(kkrsz_ns*kkrsz_ns,tmat,1,tau00,1)
  for(int i=0; i<blockSize; i++)
    for(int j=0; j<blockSize; j++)
      tau00(i,j) = tMatrices[0](i,j);
//      call zgetrs('n',kkrsz_ns,kkrsz_ns,wbig,kkrsz_ns,ipvt,tau00,
//     &           kkrsz_ns,info)
  LAPACK::zgetrs_("N", &blockSize, &blockSize, &wbig(0,0), &blockSize, ipvt, &tau00(0,0), &blockSize, &info);

// c     ----------------------------------------------------------------
// c  Redefine tau00 to be tau00-t
// c  delta is 1-t*tau00^{-1} and is calculated in gettaucl
// c  and then rotated into the local frame
// c     call scale_tau00(tau00_g,kkrsz,kkrsz,lofk,n_spin_cant,
// c    &                 kappa_rmt)

/*  
      call zgemm('n','n',kkrsz_ns,kkrsz_ns,kkrsz_ns,cone,
     &           delta,kkrsz_ns,
     >           tau00,kkrsz_ns,czero,
     &           tau00_tmp,kkrsz_ns)
*/

// c     call inv_scale_tau00(tau00_tmp,kkrsz,kkrsz,lofk,n_spin_cant,
// c    &                    kappa_rmt)
}

int main(int argc, char *argv[])
{
  int matrixType, blockSize=18, numBlocks=113;
  
  printf("Test of inversion routine for LSMS\n");
  if(argc<2)
  {
    usage(argv[0]);
    return 0;
  }
  matrixType = atoi(argv[1]);
  int n = blockSize*numBlocks;

  printf("Matrix type      : %6d\n",matrixType);
  printf("Block Size       : %6d\n",blockSize);
  printf("Number of Blocks : %6d\n",numBlocks);
  printf("Total Matrix Size: %6d\n",n);
  
  std::vector<Matrix<Complex> > tMatrices;
  Matrix<Complex> G0(n, n);
  
  tMatrices.resize(numBlocks);

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
  // add unit matrix
  // for(int i=0; i<m.n_row(); i++)
  //   m(i,i) = 1.0 + m(i,i);
  Matrix<Complex> tau00Reference(blockSize, blockSize);
  solveTau00Reference(tau00Reference, m, tMatrices, blockSize, numBlocks);

  makeType1Matrix(m, G0, tMatrices, blockSize, numBlocks);
  Matrix<Complex> tau00zblocklu(blockSize, blockSize);
  solveTau00zblocklu(tau00zblocklu, m, tMatrices, blockSize, numBlocks);

  Real d = matrixDistance(tau00Reference, tau00zblocklu);
  printf("d2 (t00Reference, tau00zblocklu) = %f\n", d);
  printf("\ntau00Reference:\n"); writeMatrix(tau00Reference);
  printf("\ntau00zblocklu:\n"); writeMatrix(tau00zblocklu);

  return 0;
}