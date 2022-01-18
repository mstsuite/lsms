module clm_mod
   use, intrinsic :: iso_fortran_env
   use, intrinsic :: iso_c_binding
   implicit none

   private

   integer, parameter :: dp = real64

   real (kind = dp), parameter :: one = 1.0d0
   real (kind = dp), parameter :: two = 2.0d0
   real (kind = dp), parameter :: zero = 0.0d0
   real (kind = dp), parameter :: four = 4.0d0
   real (kind = dp), parameter :: pi = 3.14159265358979d0
   real (kind = dp), parameter :: pi4 = FOUR * PI

   public :: calc_clm

contains

   subroutine calc_clm(lmax, clm_local) bind (C, name = "calc_clm")

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

end module clm_mod
