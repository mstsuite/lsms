c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine block_inv(a,vecs,lda,na,mp,ipvt,blk_sz,nblk,delta,
     &                     iwork,rwork,work1,alg,
     &                     idcol,iprint)
!    &                     idcol,sym_ops,iprint)
c     ================================================================
c
c     ****************************************************************
c     PURPOSE:   inverse the first block of a complex matrix: a
c                (a^{-1})_00=(a_00-delta)^{-1}
c
c     INPUT:     a,      the complex matrix to be inverted
c                vecs,   the working space
c                alg,        the algorithm of the invertion
c                            = 1, UL method
c                            = 2, LU method
c                            = 3,4,5, W.A.S. method
c
c     OUTPUT:    a,   contains the first block of the inverted
c                            matrix
c                delta
c     ****************************************************************

      implicit none
      character*3 alg_name
      integer lda,na,mp,nblk
      integer ipvt(mp),blk_sz(nblk)
      integer i,j,k,ioff
      integer alg
      integer iqmr
      integer idcol(blk_sz(1))
      integer  iprint
      integer nlim
c     parameter (nlim=80)
      parameter (nlim=80)

      real*8 time,time_direct,time_qmr
      real*8 tol
      parameter (tol=1.d-8)
c     complex*16 sym_ops(blk_sz(1),blk_sz(1),blk_sz(1)-1)
!     complex*16 sym_ops(blk_sz(1),blk_sz(1),(blk_sz(1)-1)/2)
      complex*16 a(lda,na),vecs(na-blk_sz(1),blk_sz(1)*6+6)
      complex*16 work1(blk_sz(1))
      real*8     rwork(blk_sz(1))
      integer    iwork(blk_sz(1))
      complex*16 delta(blk_sz(1),blk_sz(1))
      complex*16 czero,cone,cmone
      parameter (czero=(0.d0,0.d0))
      parameter (cone=(1.d0,0.d0))
      parameter (cmone=(-1.d0,0.d0))
c
c     ================================================================
c     Invert the KKR-Matrix using the method defined by 'alg'.........
c     If using QMR and either the method fails or it takes longer than
c     the direct method (LU or UL) then the routine uses LU thereafter
c     ================================================================
       do i=1,blk_sz(1)
       do j=1,blk_sz(1)
         a(j,i)=czero
       enddo
       enddo


      if(alg.le.2.or.alg.ge.10) then
c        =============================================================
c        Use the  LU algorithm........................................
c        =============================================================
         call zblock_lu(a,lda,blk_sz,nblk,ipvt,mp,idcol,k)
c        -------------------------------------------------------------
      else
         write(6,'('' BLOCK_INV:: incorrect alg value:'',1i5)')alg
      endif

 	do i=1,blk_sz(1)
 	do j=1,blk_sz(1)
 	  delta(j,i)=-a(j,i)
 	enddo
 	enddo
!        write(*,*) delta(1,1),delta(blk_sz(1)/2+1,blk_sz(1)/2+1)
!      if(idcol(1).eq.0) then
!	if(iprint.ge.0) then
!	  write(6,'(''block_inv:: looking for sym'')')
!	endif
!	call find_sym(delta,sym_ops,blk_sz(1),idcol,ipvt,iprint)
!      else
!	k=0
!	ioff=0
!	do i=1,blk_sz(1)
!	  if(idcol(i).eq.i) then
!	    k=k+1
!	    do j=1,blk_sz(1)
!	      delta(j,i)=-a(j,k)
!	    enddo
!	  else
!	    j=idcol(i)
!	    ioff=ioff+1
!           call zgemv('n',blk_sz(1),blk_sz(1),cone,sym_ops(1,1,ioff),
!    &           blk_sz(1),delta(1,j),1,czero,delta(1,i),1)
!	  endif
!	enddo
!      endif

      if(alg.le.2.or.alg.ge.10) then
	 alg_name = ' LU'
	 if(time.gt.0.d0) time_direct = time
      endif

      if( iprint .ge. 1 ) then
         write(6,'(''BLOCK_INV:: Using '',a3,'': Block size='',
     &     i4,'', Time='',f10.3,
     >   '' Sec'')')  alg_name,blk_sz(2),time
         call flush(6)
      endif

      return
      end
