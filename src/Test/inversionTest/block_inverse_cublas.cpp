#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"

int zblock_lu_cublas(cublasHandle_t handle, Matrix<Complex> &a, int *blk_sz, int nblk, int *ipvt, int *idcol, DeviceData &devD);

void block_inverse_cublas(cublasHandle_t handle, Matrix<Complex> &a, int *blk_sz, int nblk, Matrix<Complex> &delta, int *ipvt, int *idcol, DeviceData &devD)
{
  int k;
  
  for(int i=0; i<blk_sz[0]; i++)
    for(int j=0; j<blk_sz[0]; j++)
      a(j,i)=0.0;

  zblock_lu_cublas(handle, a, blk_sz, nblk, ipvt, idcol, devD);
  
  for(int i=0; i<blk_sz[0]; i++)
    for(int j=0; j<blk_sz[0]; j++)
      delta(j,i) = -a(j,i);
  
}

