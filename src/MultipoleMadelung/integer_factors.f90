module integer_factors_mod

   implicit none

   private


   type, public :: IntegerFactors
      integer , allocatable :: lofk(:)
      integer , allocatable :: mofk(:)
      integer , allocatable :: jofk(:)
      integer , allocatable :: lofj(:)
      integer , allocatable :: mofj(:)
      integer , allocatable :: kofj(:)
      integer , allocatable :: bofk(:)
   contains
      procedure :: init => init_integer_factors
   end type IntegerFactors


contains



   subroutine init_integer_factors(this, l0)
      class(IntegerFactors), intent(inout) :: this
      integer, intent(in) :: l0
      !
      integer  :: kmax, jmax, l, m, kl, jl, n
      !
      kmax = (l0 + 1) * (l0 + 1)
      jmax = (l0 + 1) * (l0 + 2) / 2
      !
      allocate(this % lofk(kmax))
      allocate(this % mofk(kmax))
      allocate(this % jofk(kmax))
      allocate(this % bofk(kmax))
      allocate(this % lofj(jmax))
      allocate(this % mofj(jmax))
      allocate(this % kofj(jmax))
      !  ===================================================================
      !  calculate the factors: lofk, mofk, jofk, lofj, mofj, kofj, m1m, and bofk.
      !  ===================================================================
      jl = 0
      kl = 0
      do l = 0, l0
         n = (l + 1) * (l + 2) / 2 - l
         !
         do m = -l, l, 1
            kl = kl + 1
            this % lofk(kl) = l
            this % mofk(kl) = m
            this % bofk(kl) = kl - 2 * m
            if (m>=0) then
               this % jofk(kl) = n + m
            else
               this % jofk(kl) = n - m
            endif
         enddo
         !
         do m = 0, l
            jl = jl + 1
            this % lofj(jl) = l
            this % mofj(jl) = m
            this % kofj(jl) = (l + 1) * (l + 1) - l + m
         enddo
      enddo

   end subroutine init_integer_factors




end module integer_factors_mod