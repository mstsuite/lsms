!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mdosms(zj_flag, n, dos, zz, zj, w1, pi, &
   iprint, istop)
   !     ================================================================
   !
   implicit none
   !
   character  istop*32
   character  sname*20
   parameter (sname = 'mdosms')
   !
   integer    zj_flag
   integer    n, kkrsz_loc
   integer    iprint
   integer    i
   !
   real*8     pi
   !
   complex*16 zz(n)
   complex*16 zj
   complex*16 w1(n)
   complex*16 dos
   complex*16 czero
   complex*16 ctmp
   !
   parameter (czero = (0.0d0, 0.0d0))
   !
   !     ****************************************************************
   !     calculates the density of states................................
   !     written by w.a.s.,jr sept. 1, 1992
   !     ****************************************************************
   !
   dos = czero
   kkrsz_loc = sqrt(dble(n)) + .5d0
   !     ----------------------------------------------------------------
   !  Backward sum try to cancel singularities in high l's first.
   do i = 1, n
      dos = dos + zz(i) * w1(i)
   enddo
   !       l component dos but no ss terms
   !       if(iprint.ge.100) then
   !          write(6,'(i5,2f14.8,18f9.4,'' dosd'')')
   !    >     mynod,-dos/pi,
   !    >   (-zz((i-1)*kkrsz_loc+i)*w1((i-1)*kkrsz_loc+i)/pi,i=1,min(9,n))
   !       endif
   if(zj_flag.eq.1) then
      dos = dos - zj
   endif
   dos = -dos / pi
   !
   !     ================================================================
   if(istop.eq.sname) then
      call fstop(sname)
   endif
   !
   return
end
