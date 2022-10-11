// C++ version of zblock_lu
// to be modified for use with cublas 

#include "Complex.hpp"
#include "Matrix.hpp"

extern "C" {
// void zgetrf_(int *m, int *n, Complex *a, int *lda, int *ipvt, int *info);
// void zgetrs_(const char *, int *m, int *ioff, Complex *a, int *lda, int *ipvt, Complex *b, int *ldb, int *info);
// void zgemm_(const char *, const char *, int *m, int *n, int *k, Complex *alpha, Complex *a, int *lda, Complex *b, int *ldb, Complex *beta, Complex *c, int *ldc);
}
  
// a: input matrix -> output in block 1 of a
//
// returns: k -- returns actual number of columns in the calculated inverse

int zblock_lu_cpp(Matrix<Complex> &a, int *blk_sz, int nblk, int *ipvt, int *idcol)
{
  Complex cone = 1.0;
  Complex cmone = -1.0;
  int k, lda, info;
  // total size of matrix = sum of block sizes
  lda = a.l_dim();
  int na=0;
  for(int i=0; i<nblk; i++) na+=blk_sz[i];

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
              // zcopy(na-blk_sz[0],&a(blk_sz[0],i),1,&a(blk_sz[0],k),1);
              for(int j = 0; j<na-blk_sz[0]; j++)
                a(j + blk_sz[0],k) = a(j + blk_sz[0],i);
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
          LAPACK::zgetrf_(&m, &m, &a(ioff,ioff), &lda, ipvt, &info);
          if(info!=0)
          {
            printf("zgetrf info=%d  ioff=%d\n",info,ioff);
          }
// calculate the inverse of above multiplying the row block
// blk_sz(iblk) x ioff
          LAPACK::zgetrs_("n", &m, &ioff, &a(ioff,ioff), &lda, ipvt, &a(ioff,0), &lda, &info);
          if(info!=0)
          {
            printf("zgetrs info=%d  ioff=%d\n",info,ioff);
          }
          if(iblk > 1)
          {
            int off1 = ioff-k+1;
            int off2 = na-ioff;
            BLAS::zgemm_("n", "n", &n, &off1, &off2, &cmone, &a(joff,ioff), &lda, &a(ioff,k-1), &lda, &cone, &a(joff,k-1), &lda);
            BLAS::zgemm_("n", "n", &joff, &n, &off2, &cmone, &a(0,ioff), &lda, &a(ioff,joff), &lda, &cone, &a(0,joff), &lda);
          }
        }
        int off3 = blk_sz[0]-k+1;
        int off4 = na-blk_sz[0];
        BLAS::zgemm_("n", "n", &blk_sz[0], &off3, &off4, &cmone, &a(0,blk_sz[0]), &lda, &a(blk_sz[0],k-1), &lda, &cone, &a(0,0), &lda);
      }

//     write(6,*) "k out =",k
      return blk_sz[0]-k+1;    
}

