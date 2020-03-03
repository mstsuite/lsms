// C++ version of zblock_lu
// to be modified for use with cublas 

#include <cuda_runtime.h>

#include <Complex.hpp>
#include <Matrix.hpp>

#define IDX2C(i,j,ld) (((j)*(ld))+(i))

// a: input matrix -> output in block 1 of a
//
// returns: k -- returns actual number of columns in the calculated inverse

int zblock_lu_cublas(cublasHandle_t handle, Matrix<Complex> &a, int *blk_sz, int nblk, int *ipvt, int *idcol, DeviceData &devD)
{
  int k;
  int info[1];
  int lda=a.l_dim();
  const cuDoubleComplex cone  = make_cuDoubleComplex( 1.0, 0.0);
  const cuDoubleComplex cmone = make_cuDoubleComplex(-1.0, 0.0);
  cuDoubleComplex *aAddr, *bAddr;
  cublasStatus_t cublasStat;
  // total size of matrix = sum of block sizes
  int na=0;
  for(int i=0; i<nblk; i++) na+=blk_sz[i];

  cuDoubleComplex *devA=devD.m;

  /*
#ifndef BUILDKKRMATRIX_GPU
  // copy matrix to device 
  cublasStat = cublasSetMatrix ( na, na, sizeof(cuDoubleComplex), &a(0,0), lda, devA, lda); 
#endif
  */
  
// printf("idcol[0]=%d\n",idcol[0]);
      if(idcol[0] == 0)
        k=1;
      else
      {
// eliminate columns that are equiv due to symmetry
        k=blk_sz[0];
        for(int i=blk_sz[0]-1; i>=0; i--)
        {
          if(idcol[0]==0 || idcol[i] == i+1) // i+1 due to Fortran convention
          {
            k=k-1;
            if(k!=i)
            {
              cublasZcopy(handle, na-blk_sz[0], (cuDoubleComplex*)&devA[IDX2C(blk_sz[0],i)], 1, (cuDoubleComplex*)&devA[IDX2C(blk_sz[0],k)], 1);  // check k vs k-1 ?
            }
          }
        }
      }

  

      if(nblk>0)
      {
// Do block LU
        int n=blk_sz[nblk-1];
        int joff=na-n;
        for(int iblk=nblk-1; iblk>=1; iblk--)
        {
          int m=n;
          int ioff=joff;
          n=blk_sz[iblk-1];
          joff=joff-n;
// invert the diagonal blk_sz(iblk) x blk_sz(iblk) block
          // aAddr= (cuDoubleComplex*) &a(ioff,ioff);
          aAddr= (cuDoubleComplex*) &devA[IDX2C(ioff,ioff,lda)];
          cublasZgetrfBatched(handle, m, &aAddr, lda, devD.ipvt, devD.info, 1);
          // cudaDeviceSynchronize();
          // zgetrf_(&m, &m, &a(ioff,ioff), &lda, ipvt, info); 
          if(*info!=0)
          {
            printf("zgetrf info=%d  ioff=%d\n",*info,ioff);
          }
// calculate the inverse of above multiplying the row block
// blk_sz(iblk) x ioff
          // aAddr = (cuDoubleComplex*) &a(ioff,ioff);
          // bAddr = (cuDoubleComplex*) &a(ioff,0);
          aAddr = (cuDoubleComplex*) &devA[IDX2C(ioff,ioff,lda)];
          bAddr = (cuDoubleComplex*) &devA[IDX2C(ioff,0,lda)];
          cublasZgetrsBatched(handle, CUBLAS_OP_N, m, ioff, (const cuDoubleComplex**)&aAddr, lda, devD.ipvt, &bAddr, lda, info, 1);
          // cudaDeviceSynchronize();
          // zgetrs_("n", &m, &ioff, &a(ioff,ioff), &lda, ipvt, &a(ioff,0), &lda, info);
          if(*info!=0)
          {
            printf("zgetrs info=%d  ioff=%d\n",*info,ioff);
          }
          if(iblk>1)
          {
            cublasZgemm(handle, CUBLAS_OP_N,CUBLAS_OP_N, n, ioff-k+1 ,na-ioff, &cmone, (cuDoubleComplex*)&devA[IDX2C(joff,ioff,lda)], lda,
                        (cuDoubleComplex*)&devA[IDX2C(ioff,k-1,lda)], lda, &cone, (cuDoubleComplex*)&devA[IDX2C(joff,k-1,lda)], lda);
          // cudaDeviceSynchronize();
            cublasZgemm(handle, CUBLAS_OP_N,CUBLAS_OP_N, joff, n, na-ioff, &cmone, (cuDoubleComplex*)&devA[IDX2C(0,ioff,lda)], lda,
                        (cuDoubleComplex*)&devA[IDX2C(ioff,joff,lda)], lda, &cone, (cuDoubleComplex*)&devA[IDX2C(0,joff,lda)], lda);
          // cudaDeviceSynchronize();
            // int off1 = ioff-k+1;
            // int off2 = na-ioff;
            // zgemm_("n", "n", &n, &off1, &off2, (Complex *)&cmone, &a(joff,ioff), &lda, &a(ioff,k-1), &lda, (Complex *)&cone, &a(joff,k-1), &lda);
            // zgemm_("n", "n", &joff, &n, &off2, (Complex *)&cmone, &a(0,ioff), &lda, &a(ioff,joff), &lda, (Complex *)&cone, &a(0,joff), &lda);

          }
        }
        cublasZgemm(handle, CUBLAS_OP_N,CUBLAS_OP_N, blk_sz[0], blk_sz[0]-k+1, na-blk_sz[0], &cmone, (cuDoubleComplex*)&devA[IDX2C(0,blk_sz[0],lda)], lda,
                    (cuDoubleComplex*)&devA[IDX2C(blk_sz[0],k-1,lda)], lda, &cone, (cuDoubleComplex*)&devA[IDX2C(0,0,lda)], lda);
        // cudaDeviceSynchronize();

        cublasStat = cublasGetMatrix( blk_sz[0], blk_sz[0], sizeof(cuDoubleComplex), devA, lda, &a(0,0), lda);

      }



      return k=blk_sz[0]-k+1;    
}


