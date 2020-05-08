/* -*- mode: C++; c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#include "Complex.hpp"

// we might want to distinguish between systems where all lmax (and consequently kkrsz_ns) are the same
// and systems with potential different lmax on different atoms and l steps

inline void calculateHankel(Complex prel, double r, int lend, Complex *hfn)
{
  const Complex sqrtm1(0.0, 1.0);
  Complex z=prel*r;
  hfn[0]=-sqrtm1;
  hfn[1]=-1.0-sqrtm1/z;
  for(int l=1; l<lend; l++)
  {
    hfn[l+1]=(2.0*l+1)*hfn[l]/z - hfn[l-1];
  }
  /*
c             l+1
c     hfn = -i   *h (k*R  )*sqrt(E)
c                  l    ij
*/
  z=std::exp(sqrtm1*z)/rmag;
  for(int l=0; l<=lend;l++)
  {
    hfn[l]=-hfn[l]*z*ilp1[l];     
  }
}

void buildBGijCuda((LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom, int iie,
                    Complex energy, Complex *devBgij)
{
}

void buildKKRMatrixLMaxIdenticalCuda(LSMSSystemParameters &lsms, LocalTypeInfo &local, DeviceStorage &d, AtomData &atom, int iie,
                        Complex *tMatrix, Complex *devM)
{
  cublasHandle_t cublasHandle = DeviceStorage::getCublasHandle();
  int nrmat_ns = lsms.n_spin_cant*atom.nrmat; // total size of the kkr matrix
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block

  
}
