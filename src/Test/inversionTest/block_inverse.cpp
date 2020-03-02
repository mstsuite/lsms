#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"

#if defined(ACCELERATOR_CUBLAS) || defined(ACCELERATOR_CUDA_C)
#include "Accelerator/DeviceStorage.hpp"
#endif

extern "C" {
  void zblock_lu_(Complex *a, int *lda, int *blk_sz, int *nblk, int *ipvt, int *mp, int *idcol, int *k);
  void  zblock_lu_cuda_c_(Complex *a, int *lda, int *blk_sz, int *nblk, int *ipvt, int *mp, int *idcol, int *k);
}

int zblock_lu_cpp(Matrix<Complex> &a, int *blk_sz, int nblk, int *ipvt, int *idcol);

#if defined(ACCELERATOR_CUBLAS)
int zblock_lu_cublas(cublasHandle_t handle, Matrix<Complex> &a, int *blk_sz, int nblk, int *ipvt, int *idcol);
#endif

void block_inverse(Matrix<Complex> &a, int *blk_sz, int nblk, Matrix<Complex> &delta, int *ipvt, int *idcol)
{
  int k;
  
  for(int i=0; i<blk_sz[0]; i++)
    for(int j=0; j<blk_sz[0]; j++)
      a(j,i)=0.0;

#if defined(ACCELERATOR_CUBLAS)
  zblock_lu_cublas(DeviceStorage::getCublasHandle(), a, blk_sz, nblk, ipvt, idcol);
#elif  defined(ACCELERATOR_CUDA_C)
  int mp=a.l_dim();
  int lda=a.l_dim();
  zblock_lu_cuda_c_( &a(0,0), &lda, blk_sz, &nblk, ipvt, &mp, idcol, &k);
#else
  k=zblock_lu_cpp(a, blk_sz, nblk, ipvt, idcol);
  // int mp=a.l_dim();
  // int lda=a.l_dim();
  // zblock_lu_(&a(0,0), &lda, blk_sz, &nblk, ipvt, &mp, idcol, &k);
#endif
  
  for(int i=0; i<blk_sz[0]; i++)
    for(int j=0; j<blk_sz[0]; j++)
      delta(j,i) = -a(j,i);
  
}

