/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#include "linearSolvers.hpp"

#include <stdio.h>

#include "Complex.hpp"
#include "Matrix.hpp"
#include <vector>

static void buildKKRSizeTMatrix(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie,
                                Matrix<Complex> &tMatrix, int ispin) {
  // assume Matrix<Complex> tMatrix(nrmat_ns, kkrsz_ns);
  int nrmat_ns = lsms.n_spin_cant * atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant * atom.kkrsz; // size of t00 block
  int i0 = 0; // start index of the current atom's block

  tMatrix = 0.0;

  if (lsms.n_spin_pola == lsms.n_spin_cant) { // non polarized or spin canted

    for (int i = 0; i < kkrsz_ns; i++) {
      for (int j = 0; j < kkrsz_ns; j++) {
        tMatrix(i, j) = local.tmatStore(i + j * kkrsz_ns + iie * local.blkSizeTmatStore, atom.LIZStoreIdx[0]);
      }
    }

  } else {

    int jsm = kkrsz_ns * kkrsz_ns * ispin;

    for (int i = 0; i < kkrsz_ns; i++) {
      for (int j = 0; j < kkrsz_ns; j++) {
        tMatrix(i, j) = local.tmatStore(i + j * kkrsz_ns + iie * local.blkSizeTmatStore + jsm, atom.LIZStoreIdx[0]);
      }
    }

  }

}

// given the  m-Matrix [(blockSize*numBlocks) x (blockSize*numBlocks)] and tMatrix[0] [blockSize x blockSize]
// calculate the t00-Matrix as the blockSize x blockSize diagonal block of  (1 - tG)^-1 t
void solveTau00zgesv(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie, Matrix<Complex> &m,
                     Matrix<Complex> &tau00, int ispin) {
  int nrmat_ns = lsms.n_spin_cant * atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant * atom.kkrsz; // size of t00 block

  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK
  Matrix<Complex> tau(nrmat_ns, kkrsz_ns);
  // copy t[0] into the top part of tau
  buildKKRSizeTMatrix(lsms, local, atom, iie, tau, ispin);

  int ipiv[nrmat_ns];
  int info;
  LAPACK::zgesv_(&nrmat_ns, &kkrsz_ns, &m(0, 0), &nrmat_ns, ipiv, &tau(0, 0), &nrmat_ns, &info);

  // copy result into tau00
  for (int i = 0; i < kkrsz_ns; i++)
    for (int j = 0; j < kkrsz_ns; j++)
      tau00(i, j) = tau(i, j);
}

void solveTau00zgetrf(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie, Matrix<Complex> &m,
                      Matrix<Complex> &tau00, int ispin) {
  int nrmat_ns = lsms.n_spin_cant * atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant * atom.kkrsz; // size of t00 block

  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK
  Matrix<Complex> tau(nrmat_ns, kkrsz_ns);
  // copy t[0] into the top part of tau
  buildKKRSizeTMatrix(lsms, local, atom, iie, tau, ispin);

  int ipiv[nrmat_ns];
  Matrix<Complex> work(nrmat_ns, kkrsz_ns);
  std::vector<std::complex<float> > swork(nrmat_ns * (nrmat_ns + kkrsz_ns));
  std::vector<double> rwork(nrmat_ns);
  int info, iter;

  LAPACK::zgetrf_(&nrmat_ns, &nrmat_ns, &m(0, 0), &nrmat_ns, &ipiv[0], &info);
  LAPACK::zgetrs_("N", &nrmat_ns, &kkrsz_ns, &m(0, 0), &nrmat_ns, &ipiv[0], &tau(0, 0), &nrmat_ns, &info);

  // copy result into tau00
  for (int i = 0; i < kkrsz_ns; i++)
    for (int j = 0; j < kkrsz_ns; j++)
      tau00(i, j) = tau(i, j);
}

#ifndef ARCH_IBM

void solveTau00zcgesv(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie, Matrix<Complex> &m,
                      Matrix<Complex> &tau00, int ispin) {
  int nrmat_ns = lsms.n_spin_cant * atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant * atom.kkrsz; // size of t00 block

  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK
  Matrix<Complex> tau(nrmat_ns, kkrsz_ns);
  Matrix<Complex> t(nrmat_ns, kkrsz_ns);
  // copy t[0] into the top part of t
  buildKKRSizeTMatrix(lsms, local, atom, iie, tau, ispin);
  tau = 0.0;

  int ipiv[nrmat_ns];
  Matrix<Complex> work(nrmat_ns, kkrsz_ns);
  std::vector<std::complex<float> > swork(nrmat_ns * (nrmat_ns + kkrsz_ns));
  std::vector<double> rwork(nrmat_ns);
  int info, iter;
  LAPACK::zcgesv_(&nrmat_ns, &kkrsz_ns, &m(0, 0), &nrmat_ns, ipiv, &t(0, 0), &nrmat_ns, &tau(0, 0), &nrmat_ns,
                  &work(0, 0), &swork[0], &rwork[0],
                  &iter, &info);

  // copy result into tau00
  for (int i = 0; i < kkrsz_ns; i++)
    for (int j = 0; j < kkrsz_ns; j++)
      tau00(i, j) = tau(i, j);
}

#endif

extern "C"
{
void
block_inv_(Complex *a, Complex *vecs, int *lda, int *na, int *mp, int *ipvt, int *blk_sz, int *nblk, Complex *delta,
           int *iwork, double *rwork, Complex *work1, int *alg,
           int *idcol, int *iprint);

}

void
solveTau00zblocklu_f77(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie, Matrix<Complex> &m,
                       Matrix<Complex> &tau00, int ispin) {
  int nrmat_ns = lsms.n_spin_cant * atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant * atom.kkrsz; // size of t00 block
  int ipvt[nrmat_ns];
  int info;
  int alg = 2, iprint = 0;
  Complex vecs[42]; // dummy, not used for alg=2
  int nblk = 3;
  Matrix<Complex> delta(kkrsz_ns, kkrsz_ns);
  int iwork[nrmat_ns];
  Real rwork[nrmat_ns];
  Complex work1[nrmat_ns];

  int blk_sz[1000];

  blk_sz[0] = kkrsz_ns;
  if (nblk == 1)
    blk_sz[0] = nrmat_ns;
  else if (nblk == 2)
    blk_sz[1] = nrmat_ns - blk_sz[0];
  else if (nblk > 2) {
    int min_sz = (nrmat_ns - blk_sz[0]) / (nblk - 1);
    int rem = (nrmat_ns - blk_sz[0]) % (nblk - 1);
    int i = 1;
    for (; i <= rem; i++)
      blk_sz[i] = min_sz + 1;
    for (; i < nblk; i++)
      blk_sz[i] = min_sz;
  }

  int idcol[blk_sz[0]];
  idcol[0] = 0;

  // with m = [[A B][C D]], A: blk_sz[0] x blk_sz[0]
  // calculate the Schur complement m/D of m with A set to zero,
  // i.e. delta = B D^-1 C
  block_inv_(&m(0, 0), vecs, &nrmat_ns, &nrmat_ns, &nrmat_ns, ipvt,
             blk_sz, &nblk, &delta(0, 0),
             iwork, rwork, work1, &alg, idcol, &iprint);


  Matrix<Complex> wbig(kkrsz_ns, kkrsz_ns);
// setup unit matrix...............................................
// n.b. this is the top diagonal block of the kkr matrix m
//      i.e. 1 - t_0 G_00, with G_ii == 0 this is just the unit matrix

  unitMatrix(wbig);

// c     get 1-delta and put it in wbig

  for (int i = 0; i < kkrsz_ns; i++)
    for (int j = 0; j < kkrsz_ns; j++)
      wbig(i, j) -= delta(i, j);
//  c     ================================================================
// c     create tau00 => {[1-t*G]**(-1)}*t : for central site only.......
// c     ----------------------------------------------------------------

  LAPACK::zgetrf_(&kkrsz_ns, &kkrsz_ns, &wbig(0, 0), &kkrsz_ns, ipvt, &info);

  int jsm;
  if (lsms.n_spin_pola == lsms.n_spin_cant) { // non polarized or spin canted
    jsm = 0;
  } else {
    jsm = kkrsz_ns * kkrsz_ns * ispin;
  }

  for (int i = 0; i < kkrsz_ns; i++)
    for (int j = 0; j < kkrsz_ns; j++)
      tau00(i, j) = local.tmatStore(i + j * kkrsz_ns + iie * local.blkSizeTmatStore + jsm, atom.LIZStoreIdx[0]);

  LAPACK::zgetrs_("N", &kkrsz_ns, &kkrsz_ns, &wbig(0, 0), &kkrsz_ns, ipvt, &tau00(0, 0), &kkrsz_ns, &info);
}

void block_inverse(Matrix<Complex> &a, int *blk_sz, int nblk, Matrix<Complex> &delta, int *ipvt, int *idcol);

void
solveTau00zblocklu_cpp(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie, Matrix<Complex> &m,
                       Matrix<Complex> &tau00, int ispin) {
  int nrmat_ns = lsms.n_spin_cant * atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant * atom.kkrsz; // size of t00 block
  int ipvt[nrmat_ns];
  int info;
  int alg = 2, iprint = 0;
  Complex vecs[42]; // dummy, not used for alg=2
  int nblk = 3;
  Matrix<Complex> delta(kkrsz_ns, kkrsz_ns);
  int iwork[nrmat_ns];
  Real rwork[nrmat_ns];
  Complex work1[nrmat_ns];

  int blk_sz[1000];

  blk_sz[0] = kkrsz_ns;
  if (nblk == 1)
    blk_sz[0] = nrmat_ns;
  else if (nblk == 2)
    blk_sz[1] = nrmat_ns - blk_sz[0];
  else if (nblk > 2) {
    int min_sz = (nrmat_ns - blk_sz[0]) / (nblk - 1);
    int rem = (nrmat_ns - blk_sz[0]) % (nblk - 1);
    int i = 1;
    for (; i <= rem; i++)
      blk_sz[i] = min_sz + 1;
    for (; i < nblk; i++)
      blk_sz[i] = min_sz;
  }

  int idcol[blk_sz[0]];
  idcol[0] = 0;

  // with m = [[A B][C D]], A: blk_sz[0] x blk_sz[0]
  // calculate the Schur complement m/D of m with A set to zero,
  // i.e. delta = B D^-1 C
  block_inverse(m, blk_sz, nblk, delta, ipvt, idcol);

  Matrix<Complex> wbig(kkrsz_ns, kkrsz_ns);
// setup unit matrix...............................................
// n.b. this is the top diagonal block of the kkr matrix m
//      i.e. 1 - t_0 G_00, with G_ii == 0 this is just the unit matrix

  unitMatrix(wbig);

// c     get 1-delta and put it in wbig

  for (int i = 0; i < kkrsz_ns; i++)
    for (int j = 0; j < kkrsz_ns; j++)
      wbig(i, j) -= delta(i, j);
//  c     ================================================================
// c     create tau00 => {[1-t*G]**(-1)}*t : for central site only.......
// c     ----------------------------------------------------------------

  LAPACK::zgetrf_(&kkrsz_ns, &kkrsz_ns, &wbig(0, 0), &kkrsz_ns, ipvt, &info);

  int jsm;
  if (lsms.n_spin_pola == lsms.n_spin_cant) { // non polarized or spin canted
    jsm = 0;
  } else {
    jsm = kkrsz_ns * kkrsz_ns * ispin;
  }

  for (int i = 0; i < kkrsz_ns; i++)
    for (int j = 0; j < kkrsz_ns; j++)
      tau00(i, j) = local.tmatStore(i + j * kkrsz_ns + iie * local.blkSizeTmatStore + jsm, atom.LIZStoreIdx[0]);

  LAPACK::zgetrs_("N", &kkrsz_ns, &kkrsz_ns, &wbig(0, 0), &kkrsz_ns, ipvt, &tau00(0, 0), &kkrsz_ns, &info);
}

