c
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine zeroout(x,nx)
c     =================================================================
c
      implicit   none
c
      integer    nx,n
c
      real*8     x(nx)
      real*8     zero
c
      parameter  (zero=0.0)
c
      do n=1,nx
         x(n)=zero
      enddo
c
      return
      end

c
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine czeroout(x,nx)
c     =================================================================
c
      implicit   none
c
      integer    nx,n
c
      complex*16     x(nx)
      complex*16     zero
c
      parameter  (zero=(0.0,0.0))
c
      do n=1,nx
         x(n)=zero
      enddo
c
      return
      end

