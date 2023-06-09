/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
// calculate the tau00 using the zblocklu solver

#define IDX(i, j, ldim) ((j) * (ldim) + (i))
#define IDX3(i, j, k, dim1, dim2) \
  (((k) * (dim1) * (dim2)) + ((j) * (dim1)) + (i))

// solveTau00<spin>_<solver>
// where <spin> is either SpinCanted (for non-collinear & relativistic) or
// Collinear (for spin polarized or non polarized) and <solver> specifies the
// solver used

void solveTau00SpinCanted_zblocklu(LSMSSystemParameters &lsms,
                                   LocalTypeInfo &local, AtomData &atom,
                                   Matrix<Complex> &tau00, Matrix<Complex> &m) {
  int blockSize = lsms.n_spin_cant * atom.kkrsz;  // kkrsz_ns
  int nrmat_ns = m.l_dim();                       // this should be nrmat_ns
  int ipvt[nrmat_ns];
  int info;
  int alg = 2, iprint = 0;
  Complex vecs[42];  // dummy, not used for alg=2
  int nblk = 3;
  Matrix<Complex> delta(blockSize, blockSize);
  int iwork[nrmat_ns];
  Real rwork[nrmat_ns];
  Complex work1[nrmat_ns];

  int blk_sz[1000];

  blk_sz[0] = blockSize;
  if (nblk == 1)
    blk_sz[0] = nrmat_ns;
  else if (nblk == 2)
    blk_sz[1] = nrmat_ns - blk_sz[0];
  else if (nblk > 2) {
    int min_sz = (nrmat_ns - blk_sz[0]) / (nblk - 1);
    int rem = (nrmat_ns - blk_sz[0]) % (nblk - 1);
    int i = 1;
    for (; i <= rem; i++) blk_sz[i] = min_sz + 1;
    for (; i < nblk; i++) blk_sz[i] = min_sz;
  }

  int idcol[blk_sz[0]];
  idcol[0] = 0;

  // with m = [[A B][C D]], A: blk_sz[0] x blk_sz[0]
  // calculate the Schur complement m/D of m with A set to zero,
  // i.e. delta = B D^-1 C
  block_inv_(&m(0, 0), vecs, &nrmat_ns, &nrmat_ns, &nrmat_ns, ipvt, blk_sz,
             &nblk, &delta(0, 0), iwork, rwork, work1, &alg, idcol, &iprint);

  Matrix<Complex> wbig(blockSize, blockSize);
  // setup unit matrix...............................................
  // n.b. this is the top diagonal block of the kkr matrix m
  //      i.e. 1 - t_0 G_00, with G_ii == 0 this is just the unit matrix

  unitMatrix(wbig);

  // c     get 1-delta and put it in wbig

  for (int i = 0; i < blockSize; i++)
    for (int j = 0; j < blockSize; j++) wbig(i, j) -= delta(i, j);
  //  c     ================================================================
  // c     create tau00 => {[1-t*G]**(-1)}*t : for central site only.......
  // c     ----------------------------------------------------------------

  LAPACK::zgetrf_(&blockSize, &blockSize, &wbig(0, 0), &blockSize, ipvt, &info);

  for (int i = 0; i < blockSize; i++)
    for (int j = 0; j < blockSize; j++)
      tau00(i, j) =
          local.tmatStore(iie * local.blkSizeTmatStore + IDX(i, j, blockSize),
                          atom.LIZStoreIdx[0]);  // tMatrices[0](i,j);

  LAPACK::zgetrs_("N", &blockSize, &blockSize, &wbig(0, 0), &blockSize, ipvt,
                  &tau00(0, 0), &blockSize, &info);
}

void solveTau00Collinear_zblocklu(LSMSSystemParameters &lsms,
                                  LocalTypeInfo &local, AtomData &atom,
                                  int spin, Matrix<Complex> &tau00,
                                  Matrix<Complex> &m) {
  int blockSize = lsms.n_spin_cant * atom.kkrsz;  // kkrsz_ns
  int nrmat_ns = m.l_dim();                       // this should be nrmat_ns
  int ipvt[nrmat_ns];
  int info;
  int alg = 2, iprint = 0;
  Complex vecs[42];  // dummy, not used for alg=2
  int nblk = 3;
  Matrix<Complex> delta(blockSize, blockSize);
  int iwork[nrmat_ns];
  Real rwork[nrmat_ns];
  Complex work1[nrmat_ns];

  int blk_sz[1000];

  blk_sz[0] = blockSize;
  if (nblk == 1)
    blk_sz[0] = nrmat_ns;
  else if (nblk == 2)
    blk_sz[1] = nrmat_ns - blk_sz[0];
  else if (nblk > 2) {
    int min_sz = (nrmat_ns - blk_sz[0]) / (nblk - 1);
    int rem = (nrmat_ns - blk_sz[0]) % (nblk - 1);
    int i = 1;
    for (; i <= rem; i++) blk_sz[i] = min_sz + 1;
    for (; i < nblk; i++) blk_sz[i] = min_sz;
  }

  int idcol[blk_sz[0]];
  idcol[0] = 0;

  // with m = [[A B][C D]], A: blk_sz[0] x blk_sz[0]
  // calculate the Schur complement m/D of m with A set to zero,
  // i.e. delta = B D^-1 C
  block_inv_(&m(0, 0), vecs, &nrmat_ns, &nrmat_ns, &nrmat_ns, ipvt, blk_sz,
             &nblk, &delta(0, 0), iwork, rwork, work1, &alg, idcol, &iprint);

  Matrix<Complex> wbig(blockSize, blockSize);
  // setup unit matrix...............................................
  // n.b. this is the top diagonal block of the kkr matrix m
  //      i.e. 1 - t_0 G_00, with G_ii == 0 this is just the unit matrix

  unitMatrix(wbig);

  // c     get 1-delta and put it in wbig

  for (int i = 0; i < blockSize; i++)
    for (int j = 0; j < blockSize; j++) wbig(i, j) -= delta(i, j);
  //  c     ================================================================
  // c     create tau00 => {[1-t*G]**(-1)}*t : for central site only.......
  // c     ----------------------------------------------------------------

  LAPACK::zgetrf_(&blockSize, &blockSize, &wbig(0, 0), &blockSize, ipvt, &info);

  for (int i = 0; i < blockSize; i++)
    for (int j = 0; j < blockSize; j++)
      tau00(i, j) = local.tmatStore(
          iie * local.blkSizeTmatStore +
              IDX3(i, j, spin, blockSize, blockSize * blockSize),
          atom.LIZStoreIdx[0]);  // tMatrices[0](i,j);

  LAPACK::zgetrs_("N", &blockSize, &blockSize, &wbig(0, 0), &blockSize, ipvt,
                  &tau00(0, 0), &blockSize, &info);
}
