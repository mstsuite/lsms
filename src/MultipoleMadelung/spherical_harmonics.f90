module spherical_harmonics_mod
   use, intrinsic :: iso_fortran_env
   use, intrinsic :: iso_c_binding

   implicit none

   integer, private, parameter :: dp = real64
   real (kind = dp) :: tol = 0.5 * 1.0e-12

   private

   public :: sph_harm_0
   public :: sph_harm_1


contains

   subroutine sph_harm_0(x, y, z, lmax, ylm) bind (C, name = "sph_harm_0")
      use legendre_mod, only: legendre
      use clm_mod, only: calc_clm

      real (kind = dp), intent(in) :: x, y, z
      integer, intent(in) :: lmax
      complex (kind = dp), intent(out) :: ylm((lmax + 1) * (lmax + 1))

      integer :: jmax
      integer :: kmax
      integer :: l, m, jl, kl

      complex (kind = dp), parameter :: czero = (0.0d0, 0.0d0)
      complex (kind = dp), parameter :: cone = (1.0d0, 0.0d0)
      real (kind = dp), parameter :: zero = 0.0d0
      real (kind = dp), parameter :: one = 1.0d0
      complex (kind = dp), parameter :: sqrtm1 = (0.0d0, 1.0d0)

      real(kind = dp) :: clm((lmax + 1) * (lmax + 1))
      real(kind = dp) :: plm(1:((lmax + 1) * (lmax + 2)) / 2)
      complex(kind = dp) :: e_imp(-lmax:lmax)
      integer :: m1m(-lmax : lmax)

      integer :: mofj(((lmax + 1) * (lmax + 2)) / 2)
      integer :: lofj(((lmax + 1) * (lmax + 2)) / 2)

      real (kind = dp) :: r, q, q2
      real (kind = dp) :: sin_phi
      real (kind = dp) :: cos_phi
      real (kind = dp) :: cos_the
      real (kind = dp) :: cp

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

      m1m(0)=1
      do m = 1, lmax
         m1m(m)  = -m1m(m-1)
         m1m(-m) =  m1m(m)
      enddo

      jmax = (lmax + 1) * (lmax + 2) / 2
      kmax = (lmax + 1) * (lmax + 1)
      ylm = czero
      !
      q2 = x * x + y * y
      r = sqrt(q2 + z * z)
      q = sqrt(q2)
      !
      !  ===================================================================
      !  calculating Ylm => ylm
      !  ===================================================================
      if (r<tol) then
         ylm(1) = clm(1)
         ylm(2:kmax) = czero
      else if (q<tol .and. z>zero) then   ! the point on +z-axis
         ylm(1:kmax) = czero
         do l = 0, lmax
            jl = ((l + 1) * (l + 2)) / 2 - l
            kl = (l + 1) * (l + 1) - l
            ylm(kl) = clm(jl)
         enddo
      else if (q<tol .and. z<zero) then   ! the point on -z-axis
         ylm(1:kmax) = czero
         do l = 0, lmax
            jl = ((l + 1) * (l + 2)) / 2 - l
            kl = (l + 1) * (l + 1) - l
            ylm(kl) = clm(jl) * m1m(l)
         enddo
      else
         cos_phi = x / q
         sin_phi = y / q
         cos_the = z / r
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
      !
   end subroutine sph_harm_0
   !

   subroutine sph_harm_1(vec, lmax, ylm) bind (C, name = "sph_harm_1")
      use legendre_mod, only: legendre
      use clm_mod, only: calc_clm

      real (kind = dp), intent(in) :: vec(3)
      integer, intent(in) :: lmax
      complex (kind = dp), intent(out) :: ylm((lmax + 1) * (lmax + 1))
      integer :: m1m(-lmax : lmax)

      real(kind = dp) :: clm((lmax + 1) * (lmax + 1))
      real(kind = dp) :: plm(1:((lmax + 1) * (lmax + 2)) / 2)
      complex(kind = dp) :: e_imp(-lmax:lmax)

      complex (kind = dp), parameter :: czero = (0.0d0, 0.0d0)
      complex (kind = dp), parameter :: cone = (1.0d0, 0.0d0)
      real (kind = dp), parameter :: zero = 0.0d0
      real (kind = dp), parameter :: one = 1.0d0
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

      m1m(0)=1
      do m = 1, lmax
         m1m(m)  = -m1m(m-1)
         m1m(-m) =  m1m(m)
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

end module spherical_harmonics_mod
