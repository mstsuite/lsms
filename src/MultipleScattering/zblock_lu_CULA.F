      integer function idx2f(i,j,ld)
      implicit none
      integer i,j,ld

      idx2f=((((j)-1)*(ld))+((i)-1))
      end function

      subroutine zblock_lu(a,lda,blk_sz,nblk,ipvt,mp,idcol,k)
c does a partitioning of a matrix and a partial inversion to
c get the upper subblock of the inverse
c   a -- complex*16 of (lda,size) the matrix to be inverted
c        where size is the sum of the all elements of blk_sz
c        on return, the upper subblock of a is untouched, and
c        postprocessing is needed to obtain the inverse
c   blk_sz -- integer of (nblk), the size of each subblock of matrix a
c   nblk -- number of blocks (number of elements in blk_sz)
c   ipvt -- integer of (mp), work array
c   idcol -- integer of (blk_sz(1)), if idcol(i)=idcol(j), then
c            the two columns are equivalent by symmetry, and only
c            the min(i,j) th column of the inverse is calculated.
c   k -- returns actual number of columns in the calculated inverse

C cuBLAS ZGEMM, CULA driver api GGD

      implicit none
      integer lda,na,mp,nblk
      integer ipvt(mp),blk_sz(*)
      integer i,j,k,ioff,m,joff,n,iblk,info
      integer idcol(blk_sz(1))
      complex*16 a(lda,*)
      complex*16 cone,cmone
      parameter (cone=(1.d0,0.d0))
      parameter (cmone=(-1.d0,0.d0))

!     integer sizeof_Z,sizeof_I
!     parameter (sizeof_Z=16)
!     parameter (sizeof_I=8)

      external idx2f
      integer idx2f
      external cula_device_zgetrf,cula_device_zgetrs
      external cula_zgetrf,cula_zgetrs
      external CULA_GET_ERROR_INFO,
     &         CULA_GET_ERROR_INFO_STRING,
     &         CULA_GET_LAST_STATUS
      integer cula_device_zgetrf,cula_device_zgetrs
      integer cula_zgetrf,cula_zgetrs
      integer CULA_GET_ERROR_INFO,
     &        CULA_GET_ERROR_INFO_STRING,
     &        CULA_GET_LAST_STATUS
      external cublas_alloc,cublas_set_matrix,cublas_get_matrix
      integer cublas_alloc,cublas_set_matrix,cublas_get_matrix
      integer culaStatus,culaError
      character culaErrorStr(100)
      external get_dev_m,get_dev_ipvt
      integer*8 get_dev_ipvt
      integer*8 get_dev_m
      integer*8 devA, devIPVT

      include '../Accelerator/CULA_Common.h'

      devA = get_dev_m()
      devIPVT = get_dev_ipvt()
      
      na=0
      do i=1,abs(nblk)
        na=na+blk_sz(i)
      enddo
      k=1

!      write(*,*) 'entering zblock_lu (CULA)'

!!    allocate device memory
!     info = cublas_alloc(lda*na, sizeof_Z, devA)
!     write(*,*) 'cublas_alloc of A :',info
!     info = cublas_alloc(lda,sizeof_I,devIPVT)
!     write(*,*) 'cublas_alloc of ipvt :',info
!     write(*,*) 'devA: ', devA
!     write(*,*) 'devipvt: ', devIPVT

#if !defined(BUILDKKRMATRIX_GPU)

      if(idcol(1).eq.0) then
        k=1
      else
c eliminate columns that are equiv due to symmetry
      k=blk_sz(1)+1
      do i=blk_sz(1),1,-1
        if(idcol(1).eq.0.or.idcol(i).eq.i) then
        k=k-1
        if(k.ne.i) then
        call zcopy(na-blk_sz(1),a(blk_sz(1)+1,i),1,a(blk_sz(1)+1,k),1)
        endif
        endif
      enddo
      endif

!     copy matrix to device
      info = cublas_set_matrix(lda, na, sizeof_Z, a, lda, devA, lda)
!      write(*,*) 'cublas_set_matrix of A :',info
#else
      k=1
#endif

      if(nblk.gt.0) then
c Do block LU
      n=blk_sz(nblk)
      joff=na-n
      do iblk=nblk,2,-1
!      write(*,*) 'iblk = ',iblk
      m=n
      ioff=joff
      n=blk_sz(iblk-1)
      joff=joff-n
c invert the diagonal blk_sz(iblk) x blk_sz(iblk) block
!      write(*,*) m,ioff,lda,idx2f(ioff+1,ioff+1,lda)
!     call zgetrf(m,m,a(ioff+1,ioff+1),lda,ipvt,info)
!     info = cula_zgetrf(m,m,a(ioff+1,ioff+1),lda,ipvt)
      info = cula_device_zgetrf(m,m,
     &       devA+idx2f(ioff+1,ioff+1,lda)*sizeof_Z,lda,devIPVT)
      culaStatus = cula_get_last_status()
      if (culaStatus.ne.0) then
        call CULA_GET_STATUS_STRING(culaStatus)
        culaError = CULA_GET_ERROR_INFO()
        culaStatus = CULA_GET_ERROR_INFO_STRING(culaStatus,culaError,
     &                                       culaErrorStr,100)
         write(*,*)"zgetrf error: ",culaError
         write(*,*) culaErrorStr
        STOP !JL May want more graceful shutdown
      endif
!     write(*,*) 'zgetrf info=',info
c calculate the inverse of above multiplying the row block
c blk_sz(iblk) x ioff
!     call zgetrs('n',m,ioff,a(ioff+1,ioff+1),lda,ipvt,
!    &     a(ioff+1,1),lda,info)
!     info = cula_zgetrs('n',m,ioff,a(ioff+1,ioff+1),lda,ipvt,
!    &     a(ioff+1,1),lda)
      info = cula_device_zgetrs('n',m,ioff,
     &     devA+idx2f(ioff+1,ioff+1,lda)*sizeof_Z,lda,devIPVT,
     &     devA+idx2f(ioff+1,1,lda)*sizeof_Z,lda)
      culaStatus = cula_get_last_status()
      if (culaStatus.ne.0) then
        call CULA_GET_STATUS_STRING(culaStatus)
        culaError = CULA_GET_ERROR_INFO()
        culaStatus = CULA_GET_ERROR_INFO_STRING(culaStatus,culaError,
     &                                       culaErrorStr,100)
         write(*,*)"zgetrs error: ",culaError
         write(*,*) culaErrorStr
        STOP !JL May want more graceful shutdown
      endif
      if(iblk.gt.2) then
!     call zgemm('n','n',n,ioff-k+1,na-ioff,cmone,a(joff+1,ioff+1),lda,
!    &     a(ioff+1,k),lda,cone,a(joff+1,k),lda)
!     info = cublas_set_matrix(lda, na, sizeof_Z, a, lda, devA, lda)
      call cublas_zgemm('n','n',n,ioff-k+1,na-ioff,cmone,
     &     devA+idx2f(joff+1,ioff+1,lda)*sizeof_Z,lda,
     &     devA+idx2f(ioff+1,k,lda)*sizeof_Z,lda,cone,
     &     devA+idx2f(joff+1,k,lda)*sizeof_Z,lda)
!     info = cublas_get_matrix(lda, na, sizeof_Z, devA, lda, a, lda)
!     call zgemm('n','n',joff,n,na-ioff,cmone,a(1,ioff+1),lda,
!    &     a(ioff+1,joff+1),lda,cone,a(1,joff+1),lda)
!     info = cublas_set_matrix(lda, na, sizeof_Z, a, lda, devA, lda)
      call cublas_zgemm('n','n',joff,n,na-ioff,cmone,
     &     devA+idx2f(1,ioff+1,lda)*sizeof_Z,lda,
     &     devA+idx2f(ioff+1,joff+1,lda)*sizeof_Z,lda,cone,
     &     devA+idx2f(1,joff+1,lda)*sizeof_Z,lda)
!     info = cublas_get_matrix(lda, na, sizeof_Z, devA, lda, a, lda)
      endif
      enddo
!     call zgemm('n','n',blk_sz(1),blk_sz(1)-k+1,na-blk_sz(1),cmone,
!    &     a(1,blk_sz(1)+1),lda,a(blk_sz(1)+1,k),lda,cone,a,lda)
!     info = cublas_set_matrix(lda, na, sizeof_Z, a, lda, devA, lda)
      call cublas_zgemm('n','n',blk_sz(1),blk_sz(1)-k+1,na-blk_sz(1),
     &     cmone,devA+idx2f(1,blk_sz(1)+1,lda)*sizeof_Z,lda,
     &     devA+idx2f(blk_sz(1)+1,k,lda)*sizeof_Z,lda,cone,devA,lda)
!     info = cublas_get_matrix(lda, na, sizeof_Z, devA, lda, a, lda)
!      info = cublas_get_matrix(lda,blk_sz(1),sizeof_Z,devA,lda,a,lda)
      info = 
     & cublas_get_matrix(blk_sz(1),blk_sz(1),sizeof_Z,devA,lda,a,lda)
      endif ! nblk.gt.0

!     info = cublas_get_matrix(lda, na, sizeof_Z, devA, lda, a, lda)

!     call cublas_free(devA)
!     call cublas_free(devIPVT)

!      write(*,*) 'leaving zblock_lu (CULA)'

      k=blk_sz(1)-k+1
      return
      end
