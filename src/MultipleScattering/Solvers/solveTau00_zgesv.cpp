/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
// calculate the tau00 using the full linear solver

#define IDX(i,j,ldim) ((j)*(ldim) + (i))
#define IDX3(i,j,k,dim1,dim2) (((k)*(dim1)*(dim2)) + ((j)*(dim1)) + (i))

// solveTau00<spin>_<solver>
// where <spin> is either SpinCanted (for non-collinear & relativistic) or Collinear (for spin polarized or non polarized)
// and <solver> specifies the solver used

void solveTau00SpinCanted_zgesv(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, Matrix<Complex> &tau00, Matrix<Complex> &m)
{
  int blockSize = lsms.n_spin_cant*atom.kkrsz; // kkrsz_ns
  int mSize = m.l_dim(); // this should be nrmat_ns
  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK
  Matrix<Complex> tau(mSize, blockSize);
  zeroMatrix(tau);
  // copy t[0] into the top part of tau
  for(int i=0; i<blockSize; i++)
    for(int j=0; j<blockSize; j++)
      tau(i,j) = local.tmatStore(iie*local.blkSizeTmatStore + IDX(i,j,blockSize),
                                 atom.LIZStoreIdx[0]); //tMatrices[0](i,j);

  int ipiv[mSize];
  int info;
  LAPACK::zgesv_(&mSize, &blockSize, &m(0,0), &mSize, ipiv, &tau(0,0), &mSize, &info);

  // copy result into tau00
  for(int i=0; i<blockSize; i++)
    for(int j=0; j<blockSize; j++)
      tau00(i,j) = tau(i,j);
}

void solveTau00Collinear_zgesv(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int spin,
                               Matrix<Complex> &tau00, Matrix<Complex> &m)
{
  int blockSize = atom.kkrsz; // kkrsz_ns
  int mSize = m.l_dim(); // this should be nrmat_ns
  // reference algorithm. Use LU factorization and linear solve for dense matrices in LAPACK
  Matrix<Complex> tau(mSize, blockSize);
  zeroMatrix(tau);
  // copy t[0] into the top part of tau
  for(int i=0; i<blockSize; i++)
    for(int j=0; j<blockSize; j++)
      tau(i,j) = local.tmatStore(iie*local.blkSizeTmatStore + IDX3(i,j,spin,blockSize,blockSize*blockSize),
                                 atom.LIZStoreIdx[0]); //tMatrices[0](i,j);

  int ipiv[mSize];
  int info;
  LAPACK::zgesv_(&mSize, &blockSize, &m(0,0), &mSize, ipiv, &tau(0,0), &mSize, &info);

  // copy result into tau00
  for(int i=0; i<blockSize; i++)
    for(int j=0; j<blockSize; j++)
      tau00(i,j) = tau(i,j);
}

