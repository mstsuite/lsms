
subroutine calc_clm(lmax, clm_local)
   use, intrinsic :: iso_fortran_env
   use, intrinsic :: iso_c_binding

   implicit none

   integer, parameter :: dp = real64

   real (kind = dp) :: tol = 0.5 * 1.0e-12
   real (kind = dp), parameter :: one = 1.0d0
   real (kind = dp), parameter :: two = 2.0d0
   real (kind = dp), parameter :: zero = 0.0d0
   real (kind = dp), parameter :: four = 4.0d0
   real (kind = dp), parameter :: pi = 3.14159265358979d0
   real (kind = dp), parameter :: pi4 = FOUR * PI

   integer, intent(in) :: lmax
   real (kind = dp), intent(inout) :: clm_local((lmax + 1) * (lmax + 2) / 2)

   integer :: l
   integer :: m
   integer :: m2
   integer :: i, j, sgn
   !
   real (kind = dp) :: xfac
   !
   !  ===================================================================
   !  Coefficients for complex spherical harmonics......................
   !  Calclates all the c(l,m)'s up to lmax............................
   !
   !              m       [ (2*l+1)*(l-|m|)!]
   !  c(l,m)= (-1)  * sqrt[-----------------]
   !                      [   4*pi*(l+|m|)! ]
   !
   !  ====================================================================
   if(lmax < 0) then
      print *, 'calClm', 'lmax < 0', lmax
      stop
   endif
   !
   clm_local(1) = sqrt(one / pi4)
   if (lmax < 50) then
      do l = 1, lmax
         xfac = sqrt(real((2 * l + 1), kind = dp) / pi4)
         sgn = 1
         do m = 0, l
            j = (l * (l + 1)) / 2 + m + 1
            clm_local(j) = one
            m2 = 2 * m
            !           The following code for calculating the factorial will overflow for large l and m values.
            do i = 1, m2
               clm_local(j) = (l - m + i) * clm_local(j)
            enddo
            clm_local(j) = sgn * xfac * sqrt(one / clm_local(j))
            sgn = -sgn
         enddo
      enddo
   else
      do l = 1, lmax
         xfac = sqrt(real((2 * l + 1), kind = dp) / pi4)
         sgn = 1
         do m = 0, l
            j = (l * (l + 1)) / 2 + m + 1
            m2 = 2 * m
            clm_local(j) = ZERO
            do i = 1, m2
               clm_local(j) = log(real(l - m + i, kind = dp)) + clm_local(j)
            enddo
            if (clm_local(j) > 256) then
               clm_local(j) = ZERO
            else
               clm_local(j) = sgn * xfac * exp(-clm_local(j) / TWO)
            endif
            sgn = -sgn
         enddo
      enddo
   endif
   !
end subroutine calc_clm


subroutine legendre(lmax, x, plm)
   use, intrinsic :: iso_fortran_env
   use, intrinsic :: iso_c_binding

   implicit none

   integer, parameter :: dp = real64

   real (kind = dp) :: tol = 0.5 * 1.0e-12
   real (kind = dp), parameter :: one = 1.0d0
   real (kind = dp), parameter :: two = 2.0d0
   real (kind = dp), parameter :: zero = 0.0d0
   real (kind = dp), parameter :: four = 4.0d0
   real (kind = dp), parameter :: pi = 3.14159265358979d0
   real (kind = dp), parameter :: pi4 = FOUR * PI

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

subroutine sph_harm_1(vec, lmax, ylm)
   use, intrinsic :: iso_fortran_env

   implicit none

   integer, parameter :: dp = real64

   real (kind = dp) :: tol = 0.5 * 1.0e-12
   real (kind = dp), parameter :: one = 1.0d0
   real (kind = dp), parameter :: two = 2.0d0
   real (kind = dp), parameter :: zero = 0.0d0
   real (kind = dp), parameter :: four = 4.0d0
   real (kind = dp), parameter :: pi = 3.14159265358979d0
   real (kind = dp), parameter :: pi4 = FOUR * PI

   real (kind = dp), intent(in) :: vec(3)
   integer, intent(in) :: lmax
   complex (kind = dp), intent(out) :: ylm((lmax + 1) * (lmax + 1))
   integer :: m1m(-lmax:lmax)

   real(kind = dp) :: clm((lmax + 1) * (lmax + 1))
   real(kind = dp) :: plm(1:((lmax + 1) * (lmax + 2)) / 2)
   complex(kind = dp) :: e_imp(-lmax:lmax)

   complex (kind = dp), parameter :: czero = (0.0d0, 0.0d0)
   complex (kind = dp), parameter :: cone = (1.0d0, 0.0d0)
   complex (kind = dp), parameter :: sqrtm1 = (0.0d0, 1.0d0)

   integer :: jmax
   integer :: kmax
   integer :: l, m, jl, kl

   real (kind = dp) :: r, q, q2
   real (kind = dp) :: sin_phi
   real (kind = dp) :: cos_phi
   real (kind = dp) :: cos_the
   real (kind = dp) :: cp

   integer :: mofj(((lmax + 1) * (lmax + 2)) / 2)
   integer :: lofj(((lmax + 1) * (lmax + 2)) / 2)

   complex (kind = dp) :: iphi

   call calc_clm(lmax, clm)

   jl = 0
   do l = 0, lmax
      do m = 0, l
         jl = jl + 1
         lofj(jl) = l
         mofj(jl) = m
      enddo
   end do

   m1m(0) = 1
   do m = 1, lmax
      m1m(m) = -m1m(m - 1)
      m1m(-m) = m1m(m)
   enddo

   ylm = czero
   !
   jmax = ((lmax + 1) * (lmax + 2)) / 2
   kmax = (lmax + 1) * (lmax + 1)
   !
   q2 = vec(1) * vec(1) + vec(2) * vec(2)
   r = sqrt(q2 + vec(3) * vec(3))
   q = sqrt(q2)
   !
   !  ===================================================================
   !  calculating Ylm => ylm
   !  ===================================================================
   if (r<tol) then
      ylm(1) = clm(1)
      ylm(2:kmax) = czero
   else if (q<tol .and. vec(3)>zero) then   ! the point on +z-axis
      ylm(1:kmax) = czero
      do l = 0, lmax
         jl = ((l + 1) * (l + 2)) / 2 - l
         kl = (l + 1) * (l + 1) - l
         ylm(kl) = clm(jl)
      enddo
   else if (q<tol .and. vec(3)<zero) then   ! the point on -z-axis
      ylm(1:kmax) = czero
      do l = 0, lmax
         jl = (l + 1) * (l + 2) / 2 - l
         kl = (l + 1) * (l + 1) - l
         ylm(kl) = clm(jl) * m1m(l)
      enddo
   else
      cos_phi = vec(1) / q
      sin_phi = vec(2) / q
      cos_the = vec(3) / r
      iphi = sqrtm1 * atan2(sin_phi, cos_phi)
      do m = -lmax, lmax
         e_imp(m) = exp(m * iphi)
      enddo
      e_imp(0) = cone
      !     ----------------------------------------------------------------
      call legendre(lmax, cos_the, plm)
      !     ----------------------------------------------------------------
      do jl = 1, jmax
         cp = clm(jl) * plm(jl)
         kl = (lofj(jl) + 1) * (lofj(jl) + 1) - lofj(jl)
         m = mofj(jl)
         ylm(kl + m) = cp * e_imp(m)
         ylm(kl - m) = m1m(m) * conjg(ylm(kl + m))
      enddo
   endif

end subroutine sph_harm_1
