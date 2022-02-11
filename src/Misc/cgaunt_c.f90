!
!     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine cgaunt(lmax, clm, &
   plmg, tg, wg, &
   cgnt, &
   lofk, mofk, &
   iprint, istop)
   !     =================================================================
   !
   implicit none
   !
   character  sname*20
   character  istop*28
   !
   integer    lmax
   integer    iprint
   integer    ng
   integer    n1
   integer    n2
   integer    j1, l1, m1
   integer    j2, l2, m2
   integer    l3, m3
   integer    lofk((2 * lmax + 1)**2)
   integer    mofk((2 * lmax + 1)**2)
   !
   real*8     plmg((2 * lmax + 1) * (2 * lmax + 2) / 2, 2 * lmax + 1)
   real*8     clm((2 * lmax + 1) * (2 * lmax + 2) / 2)
   real*8     tg(2 * (2 * lmax + 1))
   real*8     wg(2 * (2 * lmax + 1))
   real*8     gaunt
   real*8     cgnt(lmax + 1, (lmax + 1)**2, (lmax + 1)**2)
   real*8     ctemp
   real*8     ctol
   !
   parameter (ctol = 1.0d-14)
   parameter (sname = 'cgaunt')
   !
   !     ****************************************************************
   !     generate:
   !
   !               4pi  _  m1 _     m2* _     m3 _
   !        clll = int do Y  (o) * Y   (o) * Y  (o)
   !                0      l1       l2        l3
   !
   !                L3
   !             = C
   !                L1,L2
   !
   !     L1, L2: l = 0,1,..,lmax;  L3: l=0,1,...,2*lmax
   !     ****************************************************************
   !
   !     ================================================================
   !     generate the gaussian integration pts and weights...............
   !     ----------------------------------------------------------------
   call gauss_legendre_points(-1.d0, 1.d0, tg, wg, 2 * (2 * lmax + 1))
   !     ----------------------------------------------------------------
   !
   !     ================================================================
   !     generate the plm's at the gaussian nodes........................
   n1 = 2 * lmax + 1
   do ng = 1, n1
      ! meis: changed to normalized associated Legendre functions
      call plm_normalized(2 * lmax, tg(ng), plmg(1, ng))
   enddo
   !
   if(iprint.ge.1) then
      write(6, '('' lmax ='',1i4)')lmax
      write(6, '('' ng,tg(ng),wg(ng)'',i5,2d18.10)')&
         (ng, tg(ng), wg(ng), ng = 1, 2 * (2 * lmax + 1))
   endif
   !     ================================================================
   !     generate gaunt numbers..........................................
   !     ================================================================
   n1 = (lmax + 1)**2
   n2 = (lmax + 1)**2
   do j2 = 1, n2
      l2 = lofk(j2)
      m2 = mofk(j2)
      do j1 = 1, n1
         l1 = lofk(j1)
         m1 = mofk(j1)
         m3 = m2 - m1
         do l3 = l1 + l2, max(abs(m3), abs(l1 - l2)), -2
            ctemp = gaunt(l1, m1, l2, m2, l3, m3, 2 * lmax, 2 * lmax + 1, &
               wg, plmg, clm)
            if(abs(ctemp).lt.ctol) ctemp = 0.d0
            cgnt(l3 / 2 + 1, j1, j2) = ctemp
            !     ================================================================
            if(iprint.ge.2) then
               write(6, '('' l1,m1,l2,m2,l3,m3,clll = '',&
                  6i3,2x,1pd20.13)')&
                  l1, m1, &
                  l2, m2, &
                  l3, m3, &
                  cgnt(l3 / 2 + 1, j1, j2)
            endif
            !     ================================================================
            !
         enddo
      enddo
   enddo
   !
   !      if(istop.eq.sname) then
   !         call fstop(sname)
   !      else
   !         return
   !      endif
   !
end
