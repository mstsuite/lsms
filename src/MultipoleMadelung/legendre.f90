module legendre_mod
   use, intrinsic :: iso_fortran_env
   use, intrinsic :: iso_c_binding
   implicit none
   integer, parameter :: dp = real64

contains
   subroutine legendre(lmax, x, plm) bind (C, name = "legendre")

      integer, intent(in) :: lmax
      real (kind = dp), intent(in) :: x
      real (kind = dp), intent(inout) :: plm((lmax + 1) * (lmax + 2) / 2)

      integer :: jmax, l, m, i
      real (kind = dp) :: pmm
      real (kind = dp) :: somx2
      real (kind = dp) :: fact
      !
      !  ===================================================================
      !  Calclates Associated Legendre function, p(l,m), up to lmax.
      !  Based on the formulae given in "Numerical Recipes" pages 180-183
      !  (Equations 6.6.7, 6.6.8 and 6.6.9)
      !  W. H. Press, B. P. Flannery, S A Teukolsky and W. T. Vetterling.
      !  Cambridge Univ Press 1986.
      !
      !  N.B. The definition of p(l,m) has been modified.
      !  p(l,m) of this code = [(-1)**m]*p(l,m) of "Numerical Recipes".
      !  ===================================================================
      if(lmax<0) then
         print *, "Legendre", "bad lmax: lmax < 0", lmax
         stop
      else if(abs(x)>1.0_dp) then
         print *, "Legendre", "bad arguments: abs(x) > 1", x
         stop
      endif
      !
      plm = 0.0_dp
      jmax = (lmax + 1) * (lmax + 2) / 2
      if((1.0_dp - abs(x))<=1.0d-08) then
         plm(1:jmax) = 0.0_dp
         if(x<0.0_dp) then
            do l = 0, lmax
               i = (l + 1) * (l + 2) / 2 - l
               plm(i) = 1.0_dp - 2.0_dp * mod(l, 2)
            enddo
         else
            do l = 0, lmax
               i = (l + 1) * (l + 2) / 2 - l
               plm(i) = 1.0_dp
            enddo
         endif
         return
      endif
      !
      !  ===================================================================
      !  begin calculation of p(l,m)'s...................................
      !  ===================================================================
      if(lmax.eq.0) then
         !     ================================================================
         !     special case lmax=0..........................................
         !     ================================================================
         plm(1) = 1.0_dp
      else if(lmax.eq.1) then
         !     ================================================================
         !     special case lmax=1..........................................
         !     ================================================================
         plm(1) = 1.0_dp
         plm(2) = x
         plm(3) = sqrt((1.0_dp - x) * (1.0_dp + x))
      else
         plm(1) = 1.0_dp
         plm(2) = x
         somx2 = sqrt((1.0_dp - x) * (1.0_dp + x))
         do m = 1, lmax - 1
            !        =============================================================
            !                                 m       m
            !        calculate the first 2.0_dp P   and P
            !                                 m       m+1
            !        =============================================================
            pmm = 1.0_dp
            fact = 1.0_dp
            do i = 1, m
               pmm = pmm * fact * somx2
               fact = fact + 2.0_dp
            enddo
            plm(m * (m + 1) / 2 + m + 1) = pmm
            plm((m + 1) * (m + 2) / 2 + m + 1) = x * (2 * m + 1) * pmm
         enddo
         pmm = 1.0_dp
         fact = 1.0_dp
         do i = 1, lmax
            pmm = pmm * fact * somx2
            fact = fact + 2.0_dp
         enddo
         plm(lmax * (lmax + 1) / 2 + lmax + 1) = pmm
         !     ================================================================
         !                         m        m
         !     calculate the rest P     to P
         !                         m+2      lmax
         !     ================================================================
         do m = 0, lmax
            do l = m + 2, lmax
               plm(l * (l + 1) / 2 + m + 1) = (x * (2 * l - 1) * plm((l - 1) * l / 2 + m + 1) - &
                  (l + m - 1) * plm((l - 2) * (l - 1) / 2 + m + 1)) / dble(l - m)
            enddo
         enddo
      endif
      !
   end subroutine legendre

end module legendre_mod
