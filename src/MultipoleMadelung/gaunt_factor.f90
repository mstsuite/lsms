module gaunt_factor_mod
   use, intrinsic :: iso_fortran_env
   use, intrinsic :: iso_c_binding
   implicit none

   ! Parameters
   integer, parameter :: dp = real64

   real (kind = dp), parameter :: tol = 1e-14_dp
   real (kind = dp), parameter :: one = 1.0d0
   real (kind = dp), parameter :: two = 2.0d0
   real (kind = dp), parameter :: zero = 0.0d0
   real (kind = dp), parameter :: four = 4.0d0
   real (kind = dp), parameter :: ten2m14 = 1.0d-014
   real (kind = dp), parameter :: pi = 3.14159265358979d0
   real (kind = dp), parameter :: pi4 = FOUR * PI


contains

   subroutine gaunt_factor(lmax, cgnt, kj3, nj3) bind (C, name = "gaunt_factor")
      use clm_mod, only: calc_clm
      use integer_factors_mod, only: IntegerFactors

      integer, intent(in) :: lmax

      real (kind = dp), intent(inout) :: cgnt(lmax + 1, (lmax + 1)**2, (lmax + 1)**2)

      integer, intent(inout) :: kj3(lmax + 1, (lmax + 1)**2, (lmax + 1)**2)

      integer, intent(inout) :: nj3((lmax + 1)**2, (lmax + 1)**2)

      ! Auxiliar variables
      integer :: jmax, kmax, lmax2, jmax2, kmax2

      ! Common variables
      integer :: MaxJ3

      real (kind = dp) :: tol = ten2m14
      real (kind = dp), allocatable, target :: clm(:)

      type(IntegerFactors) :: integer_factors

      ! Init part starts
      kmax = (lmax + 1)**2
      jmax = (lmax + 1) * (lmax + 2) / 2
      !
      MaxJ3 = lmax + 1
      !
      lmax2 = 2 * lmax
      kmax2 = (lmax2 + 1)**2
      jmax2 = (lmax2 + 1) * (lmax + 1)
      !
      !  -------------------------------------------------------------------
      nj3(1:kmax, 1:kmax) = 0
      kj3(1:MaxJ3, 1:kmax, 1:kmax) = 0
      cgnt(1:MaxJ3, 1:kmax, 1:kmax) = 0.0
      !  -------------------------------------------------------------------
      !  set up factors lofk, mofk, etc..., for l <= lmax2
      !  -------------------------------------------------------------------
      call integer_factors % init(lmax2)
      !  -------------------------------------------------------------------
      !  get prefators of complex spherical harmonics for l <= lmax2
      !  -------------------------------------------------------------------
      allocate(clm((lmax2 + 1) * (lmax2 + 2) / 2))

      call calc_clm(lmax2, clm)

      !  -------------------------------------------------------------------
      !  set up gaunt numbers.[N.B. this routine needs the clm above],
      !
      !                   L       4pi  ->   ->   * ->     ->
      !  cgnt(j,L',L") = C      = int do*Y (o )*Y (o )*Y (o )
      !                   L',L"    0      L'     L"     L
      !
      !  L', L" <= lmax
      !
      !  L = kj3(j,L',L"), j=1,2,...,nj3(L',L") <= lmax+1
      !  -------------------------------------------------------------------

      call calc_gaunt_factor(lmax, lmax, maxj3, nj3, &
         integer_factors % lofk, &
         integer_factors % mofk, clm, cgnt, kj3)

   end subroutine gaunt_factor


   subroutine calc_gaunt_factor(lmax_1, lmax_2, MaxJ3, nj3, lofk, mofk, clm, cgnt, kj3)
      use legendre_mod, only: legendre
      use gaussq_mod, only: gaussq
      !
      implicit none
      !
      character (len = 11), parameter :: sname = 'calc_gaunt_factor'
      !
      integer, intent(in) :: lmax_1
      integer, intent(in) :: lmax_2
      integer, intent(in) :: maxj3

      integer, intent(inout) :: nj3(:, :)
      integer, intent(in) :: lofk(:)
      integer, intent(in) :: mofk(:)
      integer, intent(inout) :: kj3(:, :, :)

      real (kind = dp), intent(inout) :: clm(:)
      real (kind = dp), intent(inout) :: cgnt(:, :, :)

      integer :: lmax_3
      integer :: kmax_1
      integer :: kmax_2
      integer :: jmax_3
      integer :: kmax_3
      integer :: lm3p1
      integer :: ng
      integer :: h3
      integer :: kl1, l1, m1
      integer :: kl2, l2, m2
      integer :: l3, m3
      integer :: j3count
      !
      real (kind = dp) :: endpts(2)
      real (kind = dp) :: ctemp
      !
      real (kind = dp), allocatable :: plmg(:, :)
      real (kind = dp), allocatable :: tg(:)
      real (kind = dp), allocatable :: wg(:)
      !
      data       endpts/-1.d0, +1.d0/
      !
      lmax_3 = lmax_1 + lmax_2
      kmax_1 = (lmax_1 + 1)**2
      kmax_2 = (lmax_2 + 1)**2
      kmax_3 = (lmax_3 + 1)**2
      jmax_3 = (lmax_3 + 1) * (lmax_3 + 2) / 2
      !
      if (min(lmax_1, lmax_2) + 1 > MaxJ3) then
         print *, sname, '1st dim of cgnt is out of bound'
      endif
      !
      lm3p1 = lmax_3 + 1
      allocate(plmg(jmax_3, lm3p1))
      allocate(tg(2 * lm3p1))
      allocate(wg(2 * lm3p1))
      !
      !  ===================================================================
      !  generate the gaussian integration pts and weights...............
      !  -------------------------------------------------------------------
      call gaussq(1, 2 * lm3p1, 0, endpts(1), endpts(2), tg, wg)
      !  -------------------------------------------------------------------
      !
      !  ===================================================================
      !  generate the plm's at the gaussian nodes........................
      !  ===================================================================
      do ng = 1, lm3p1
         !    -----------------------------------------------------------------
         call legendre(lmax_3, tg(ng), plmg(1:jmax_3, ng))
         !    -----------------------------------------------------------------
      enddo
      !
      !  ===================================================================
      !
      !         4pi  _  m1 _     m2* _     m3 _
      !  clll = int do Y  (o) * Y   (o) * Y  (o)
      !          0      l1       l2        l3
      !
      !          L3
      !       = C
      !          L1,L2
      !
      !      l1 = 0, 1, ..., lmax_1
      !      l2 = 0, 1, ..., lmax_2
      !      l3 = 0, 1, ..., lmax_1+lmax_2
      !
      !  store clll in cgnt(j,L1,L2), j=1,2,...,nj3(L1,L2), L3=kj3(j,L1,L2)
      !  ===================================================================
      do kl2 = 1, kmax_2
         l2 = lofk(kl2)
         m2 = mofk(kl2)
         do kl1 = 1, kmax_1
            l1 = lofk(kl1)
            m1 = mofk(kl1)
            m3 = m2 - m1
            h3 = max(abs(m3), abs(l1 - l2))
            j3count = 0
            do l3 = l1 + l2, h3, -2
               !           ----------------------------------------------------------
               ctemp = clll(l1, m1, l2, m2, l3, m3, wg, lm3p1, plmg, lmax_3, clm)
               !           ----------------------------------------------------------
               if(abs(ctemp)>tol) then
                  if(j3count < MaxJ3) then
                     j3count = j3count + 1
                     cgnt(j3count, kl1, kl2) = ctemp
                     kj3(j3count, kl1, kl2) = (l3 + 1) * (l3 + 1) - l3 + m3
                  else
                     print *, 'j3count > MaxJ3', j3count + 1, MaxJ3
                  endif
               endif
            enddo
            nj3(kl1, kl2) = j3count
            !        do j3count=1,nj3(kl1,kl2)
            !           kl3=kj3(j3count,kl1,kl2)
            !           write(6,'('' l1,m1,l2,m2,j3,cgnt ='',5i5,1d15.8)') &
            !                        lofk(kl1),mofk(kl1),lofk(kl2),mofk(kl2), &
            !                        j3count,cgnt(j3count,kl1,kl2)
            !        enddo
         enddo
      enddo
      !
      deallocate(plmg)
      deallocate(tg)
      deallocate(wg)
      !
   end subroutine calc_gaunt_factor


   function clll(l1, m1, l2, m2, l3, m3, wg, ngauss, plmg, lmax, clm)

      integer, intent(in) :: l1
      integer, intent(in) :: l2
      integer, intent(in) :: l3
      integer, intent(in) :: m1
      integer, intent(in) :: m2
      integer, intent(in) :: m3
      integer, intent(in) :: lmax
      integer, intent(in) :: ngauss
      real (kind = dp), intent(in) :: wg(ngauss)
      real (kind = dp), intent(in) :: plmg((lmax + 1) * (lmax + 2) / 2, ngauss)
      real (kind = dp), intent(inout) :: clm(:)

      integer :: ng
      integer :: ifac1
      integer :: jl1
      integer :: ifac2
      integer :: jl2
      integer :: ifac3
      integer :: jl3
      integer :: m
      integer, target, allocatable :: m1m(:)

      real (kind = dp) :: clll

      !
      !  *******************************************************************
      !
      !               4pi  _  m1 _     m2* _     m3 _
      !        clll = int do Y  (o) * Y   (o) * Y  (o)
      !                0      l1       l2        l3
      !
      !                L3
      !             = C
      !                L1,L2
      !
      !  *******************************************************************
      !

      allocate(m1m(-lmax:lmax))

      m1m(0) = 1
      do m = 1, lmax
         m1m(m) = -m1m(m - 1)
         m1m(-m) = m1m(m)
      enddo

      if(l1.gt.lmax .or. l2.gt.lmax .or. l3.gt.lmax) then
         print *, 'CLLL', 'bad parameters: l1,l2,l3,lmax', l1, l2, l3, lmax
         stop
      endif
      clll = zero
      if(mod(l1 + l2 + l3, 2).ne.0) then
         return
      else if(m1 - m2 + m3 .ne. 0) then
         return
      else if(l1 + l2.lt.l3 .or. l2 + l3.lt.l1 .or. l3 + l1.lt.l2) then
         return
      else
         !     ----------------------------------------------------------------
         call prefac(l1, m1, ifac1, jl1, m1m, lmax)
         call prefac(l2, -m2, ifac2, jl2, m1m, lmax)
         call prefac(l3, m3, ifac3, jl3, m1m, lmax)
         !     ----------------------------------------------------------------
         do ng = 1, ngauss
            clll = clll + wg(ng) * plmg(jl1, ng) * plmg(jl2, ng) * plmg(jl3, ng)
         enddo
         if (abs(clll).lt.tol) then
            clll = zero
            return
         else
            clll = pi4 * ifac1 * clm(jl1) * m1m(m2) * ifac2 * clm(jl2) * ifac3 * clm(jl3) * clll
         endif
      endif

   end function clll


   subroutine prefac(l, m, ifac, jl, m1m, lmax)

      integer, intent(in) :: l
      integer, intent(in) :: m
      integer, intent(out) :: ifac
      integer, intent(out) :: jl
      integer, intent(in) :: lmax
      integer, intent(inout) :: m1m(-lmax:lmax)

      if (m>=0) then
         ifac = 1
         jl = (l + 1) * (l + 2) / 2 - l + m
      else
         ifac = m1m(m)
         jl = (l + 1) * (l + 2) / 2 - l - m
      end if

   end subroutine prefac

end module gaunt_factor_mod

