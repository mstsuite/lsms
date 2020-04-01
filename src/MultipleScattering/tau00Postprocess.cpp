/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */

#include "Complex.hpp"
#include "Matrix.hpp"

/*
subroutine tau_inv_postproc_nrel(kkrsz_ns,n_spin_cant,
     &                            wbig,delta,tmat,ipvt,tau00,
     &                            ubr,ubrd,
     &                            tau00_l)
      implicit none

      integer n_spin_cant
      integer kkrsz_ns
      integer ipvt(*)
      complex*16 wbig(*)
      complex*16 tmat(*)
      complex*16 tau00(*)
      complex*16 tau00_tmp(kkrsz_ns,kkrsz_ns)
      complex*16 tau00_l(*)
      complex*16 delta(*),ubr(*),ubrd(*)
      integer kkrsz

      integer info
      integer mtxsize

      complex*16 cmone,cone,czero
      parameter (cmone=(-1.0d0,0.0d0))
      parameter (cone=(1.0d0,0.0d0))
      parameter (czero=(0.d0,0.d0))

      kkrsz=kkrsz_ns/2
      mtxsize=kkrsz_ns*kkrsz_ns
!
!     FROM LSMS_1.9: GETTAU_CL
!
c     setup unit matrix...............................................
c     ----------------------------------------------------------------
++Q      call cmtruni(wbig,kkrsz_ns)
c     ----------------------------------------------------------------
c     get 1-delta and put it in wbig
++Q      call zaxpy(mtxsize,cmone,delta,1,wbig,1)
c     ================================================================
c     create tau00 => {[1-t*G]**(-1)}*t : for central site only.......
c     ----------------------------------------------------------------
++Q      call zgetrf(kkrsz_ns,kkrsz_ns,wbig,kkrsz_ns,ipvt,info)
++Q      call zcopy(kkrsz_ns*kkrsz_ns,tmat,1,tau00,1)
++Q      call zgetrs('n',kkrsz_ns,kkrsz_ns,wbig,kkrsz_ns,ipvt,tau00,
++Q     &           kkrsz_ns,info)
*/

void calculateTau00MinusT(LSMSSystemParameters &lsms, LocalTypeInfo &local, AtomData &atom, int iie, Matrix<Complex> &tau00, Matrix<Complex> &tau00MinusT)
{
  // Redefine tau00 to be tau00-t
  // delta is 1-t*tau00^{-1} and is calculated in gettaucl
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz; // size of t00 block

  //ensure that the dimensions are identical
  for(int i=0; i<kkrsz_ns; i++)
    for(int j=0; j<kkrsz_ns; j++)
      tau00MinusT(i,j) = tau00(i,j) - local.tmatStore(i + j*kkrsz_ns + iie*local.blkSizeTmatStore,atom.LIZStoreIdx[0]);
}

extern "C"
{
  void trgtol_(int *kkrsz0, int *kkrsz1, Complex *ubr, Complex *ubrd,Complex *tau00_tmp,Complex *tau00_l);
}

void rotateTau00ToLocalFrameNonRelativistic(LSMSSystemParameters &lsms, AtomData &atom, Matrix<Complex> &tau00,  Complex *tau00_l)
{
  if (lsms.n_spin_cant==2)
  {
    trgtol_(&atom.kkrsz, &atom.kkrsz, atom.ubr, atom.ubrd, &tau00(0,0), tau00_l);
  } else {
    for(int i=0; i<tau00.l_dim()*tau00.n_col(); i++)
      tau00_l[i] = tau00[i];
  }
}

/*
  // and then rotated into the local frame

c
c     ================================================================
c     Rotate tau00 to local frame of reference
c     ================================================================
      if( n_spin_cant .eq. 2 ) then
!        Non relativistic
c        -------------------------------------------------------------
         call trgtol(kkrsz,kkrsz,ubr,ubrd,tau00_tmp,tau00_l)
      else
c        -------------------------------------------------------------
         call zcopy(kkrsz_ns*kkrsz_ns,tau00_tmp,1,tau00_l,1)
c        -------------------------------------------------------------
      endif


!      write(*,*) 'tau00_l(1,1)=',tau00_l(1)

      end subroutine
*/

/*
      subroutine tau_inv_postproc_rel(kkrsz_ns,
     &                            wbig,delta,tmat,ipvt,tau00,
     &                            dmatp,dmat,
     &                            tau00_l)
      implicit none

      integer kkrsz_ns
      integer ipvt(*)
      complex*16 wbig(*)
      complex*16 tmat(*)
      complex*16 tau00(*)
      complex*16 tau00_l(*)
      complex*16 delta(*)
      complex*16 dmatp,dmat

      integer info
      integer mtxsize

      complex*16 cmone,cone,czero
      parameter (cmone=(-1.0d0,0.0d0))
      parameter (cone=(1.0d0,0.0d0))
      parameter (czero=(0.d0,0.d0))

      mtxsize=kkrsz_ns*kkrsz_ns
!
!     FROM LSMS_1.9: GETTAU_CL
!
c     setup unit matrix...............................................
c     ----------------------------------------------------------------
++Q      call cmtruni(wbig,kkrsz_ns)
c     ----------------------------------------------------------------
c     get 1-delta and put it in wbig
++Q      call zaxpy(mtxsize,cmone,delta,1,wbig,1)
c     ================================================================
c     create tau00 => {[1-t*G]**(-1)}*t : for central site only.......
c     ----------------------------------------------------------------
++Q      call zgetrf(kkrsz_ns,kkrsz_ns,wbig,kkrsz_ns,ipvt,info)
++Q      call zcopy(kkrsz_ns*kkrsz_ns,tmat,1,tau00,1)
++Q      call zgetrs('n',kkrsz_ns,kkrsz_ns,wbig,kkrsz_ns,ipvt,tau00,
++Q     &           kkrsz_ns,info)
c     ----------------------------------------------------------------

c     ================================================================
c     Rotate tau00 to local frame of reference

!        Relativistic
c        -------------------------------------------------------------
!      call zcopy(4*kkrsz_ns*kkrsz_ns,tau00,1,tau00_l,1)
      call zcopy(kkrsz_ns*kkrsz_ns,tau00,1,tau00_l,1)
      call tripmt(dmatp,tau00_l,dmat,kkrsz_ns,kkrsz_ns,kkrsz_ns)

      end subroutine
*/
extern "C"
{
  void tripmt_(Complex *dmatp,Complex *tau00_l,Complex *dmat,int *kkrsz_ns, int *, int *);
}
  
void rotateTau00ToLocalFrameRelativistic(LSMSSystemParameters &lsms, AtomData &atom, Matrix<Complex> &tau00,  Complex *tau00_l)
{
  int kkrsz_ns = lsms.n_spin_cant*atom.kkrsz;
  for(int i=0; i<tau00.l_dim()*tau00.n_col(); i++)
      tau00_l[i] = tau00[i];
  tripmt_(&atom.dmatp(0,0), tau00_l, &atom.dmat(0,0), &kkrsz_ns, &kkrsz_ns, &kkrsz_ns);
}
