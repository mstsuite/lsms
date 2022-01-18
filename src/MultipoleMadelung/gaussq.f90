module gaussq_mod
   use, intrinsic :: iso_fortran_env
   use, intrinsic :: iso_c_binding
   implicit none

   integer, parameter :: dp = real64

   real (kind = dp), parameter :: pi = 3.14159265358979d0

   private

   public :: gaussq

contains

   subroutine gaussq(kind, n, kpts, x1, x2, t, w)  bind (C, name = "gaussq")
      !  ===================================================================
      !
      !        this set of routines computes the nodes t(j) and weights
      !        w(j) for gaussian-type quadrature rules with pre-assigned
      !        nodes (x1, x2).  these are used when 1.0_dp wishes to approximate
      !
      !                 integral (from a to b)  f(x) w(x) dx
      !
      !                              n
      !        by                   sum w  f(t )
      !                             j=1  j    j
      !
      !        (note w(x) and w(j) have no connection with each other.)
      !        here w(x) is 1.0_dp of 4.0_dp possible non-negative weight
      !        functions (listed below), and f(x) is the
      !        function to be integrated.  gaussian quadrature is particularly
      !        useful on infinite intervals (with appropriate weight
      !        functions), since then other techniques often fail.
      !
      !           associated with each weight function w(x) is a set of
      !        orthogonal polynomials.  the nodes t(j) are just the 0.0_dpes
      !        of the proper n-th degree polynomial.
      !
      !        input parameters (all real numbers are in double precision)
      !
      !        kind     an integer between 1 and 6 giving the type of
      !                 quadrature rule:
      !
      !        kind = 1:  legendre quadrature, w(x) = 1 on (-1, 1)
      !        kind = 2:  chebyshev quadrature of the first kind
      !                   w(x) = 1/sqrt(1 - x*x) on (-1, +1)
      !        kind = 3:  chebyshev quadrature of the second kind
      !                   w(x) = sqrt(1 - x*x) on (-1, 1)
      !        kind = 4:  hermite quadrature, w(x) = exp(-x*x) on
      !                   (-infinity, +infinity)
      !        n        the number of points used for the quadrature rule
      !        kpts     (integer) normally 0, unless the left or right end-
      !                 point (or both) of the interval is required to be a
      !                 node (this is called gauss-radau or gauss-lobatto
      !                 quadrature).  then kpts is the number of fixed
      !                 endpoints (1 or 2).
      !        endpts   real array of length 2.  contains the values of
      !                 any fixed endpoints, if kpts = 1 or 2.
      !                 Note: now, x1 = endpts(1) and x2 = endpts(2)
      !        b        real scratch array of length n
      !                 Note: now, b is allocated internally
      !
      !        output parameters (both double precision arrays of length n)
      !
      !        t        will contain the desired nodes.
      !        w        will contain the desired weights w(j).
      !
      !        underflow may sometimes occur, but is harmless.
      !
      !        references
      !        1.  golub, g. h., and welsch, j. h., "calculation of gaussian
      !            quadrature rules," mathematics of computation 23 (april,
      !            1969), pp. 221-230.
      !        2.  golub, g. h., "some modified matrix eigenvalue problems,"
      !            siam review 15 (april, 1973), pp. 318-334 (section 7).
      !        3.  stroud and secrest, gaussian quadrature formulas, prentice-
      !            hall, englewood cliffs, n.j., 1966.
      !
      !        original version 20 jan 1975 from stanford
      !        modified 21 dec 1983 by eric grosse
      !          imtql2 => gausq2
      !          compute pi using atan
      !          removed accuracy claims, description of method
      !          added single precision version
      !        modified 15 Sept. 1991 by Yang Wang:
      !          use generic name for functions rather than the data type
      !          specified name.
      !          for simplicity, delete the routines for kind = 5 and 6,
      !          which give the roots and weights for Laguerre polynormials
      !          and Jacobi polynormials.
      !        modified 25 June 2002 by Yang Wang:
      !          b array is allocated and deallocated internally
      !          use x1 and x2 instead of endpts(1) and endpts(2)
      !  *******************************************************************
      !
      integer, intent(in) :: kind
      integer, intent(in) :: n
      integer, intent(in) :: kpts
      integer :: i, ierr
      !
      !ywg  real (kind=dp), intent(in) :: endpts(2)
      real (kind = dp), intent(in) :: x1, x2
      real (kind = dp) :: endpts(2)
      !ywg  real (kind=dp), intent(inout) :: b(n)
      real (kind = dp), allocatable :: b(:)
      real (kind = dp), intent(inout) :: t(n)
      real (kind = dp), intent(inout) :: w(n)
      real (kind = dp) :: muzero, t1, gam, xm, xl
      !
      if (n < 2) then
         print *, 'gaussq', 'n < 2', n
      endif
      allocate(b(n))
      endpts(1) = x1
      endpts(2) = x2
      !
      !  -------------------------------------------------------------------
      call class (kind, n, b, t, muzero)
      !  -------------------------------------------------------------------
      !
      !  The matrix of coefficients is assumed to be symmetric. The array t
      !  contains the diagonal elements, the array b the off-diagonal elements.
      !  Make appropriate changes in the lower right 2 by 2 submatrix.
      !
      if (kpts.eq.1) then
         !     ================================================================
         !     only t(n) must be changed
         !     ----------------------------------------------------------------
         t(n) = solve(endpts(1), n, t, b) * b(n - 1)**2 + endpts(1)
         !     ----------------------------------------------------------------
      else if (kpts.eq.2) then
         !     ================================================================
         !     t(n) and b(n-1) must be recomputed
         !     ----------------------------------------------------------------
         gam = solve(endpts(1), n, t, b)
         t1 = ((endpts(1) - endpts(2)) / (solve(endpts(2), n, t, b) - gam))
         !     ----------------------------------------------------------------
         b(n - 1) = sqrt(t1)
         t(n) = endpts(1) + gam * t1
      else if (kpts.ne.0) then
         print *, 'GAUSSQ', 'kpts < 0, or k > 2', kpts
      endif
      !
      !  note that the indices of the elements of b run from 1 to n-1
      !  and thus the value of b(n) is arbitrary.
      !  now compute the eigenvalues of the symmetric tridiagonal
      !  matrix, which has been modified as necessary.
      !  the method used is a ql-type method with origin shifting
      !
      w(1) = 1.0_dp
      w(2:n) = 0.0_dp
      !
      !  -------------------------------------------------------------------
      call gausq2 (n, t, b, w, ierr)
      !  -------------------------------------------------------------------
      if (ierr.ne.0) then
         print *, 'gaussq', 'Choose a proper machine parameter.'
         stop
      endif

      !
      if (kind /= 4) then
         xm = 0.5_dp * (x2 + x1)
         xl = 0.5_dp * (x2 - x1)
         do i = 1, n
            t(i) = xm + xl * t(i)
         enddo
         do i = 1, n
            w(i) = muzero * w(i) * w(i) * xl
         enddo
      endif



      !
      deallocate(b)
      !
   end subroutine gaussq
   !  ===================================================================
   !
   !  *******************************************************************
   !

   function solve(shift, n, a, b)
      !  ===================================================================
      !
      !       this procedure performs elimination to solve for the
      !       n-th component of the solution delta to the equation
      !
      !             (jn - shift*identity) * delta  = en,
      !
      !       where en is the vector of all zeroes except for 1 in
      !       the n-th position.
      !
      !       the matrix jn is symmetric tridiagonal, with diagonal
      !       elements a(i), off-diagonal elements b(i).  this equation
      !       must be solved to obtain the appropriate changes in the lower
      !       2 by 2 submatrix of coefficients for orthogonal polynomials.
      !  *******************************************************************

      !
      !
      !
      integer, intent(in) :: n
      integer :: i, nm1
      !
      real (kind = dp), intent(in) :: shift
      real (kind = dp), intent(in) :: a(n)
      real (kind = dp), intent(in) :: b(n)
      real (kind = dp) :: solve, alpha
      !
      alpha = a(1) - shift
      if (alpha.lt.1.0e-8_dp) then
         print *, 'SOLVE', 'Incorrect endpts or kpts'
      endif
      nm1 = n - 1
      do i = 2, nm1
         alpha = a(i) - shift - b(i - 1) * b(i - 1) / alpha
      enddo
      solve = 1.0_dp / alpha

   end function solve
   !
   !

   subroutine class(kind, n, b, a, muzero)
      !  ===================================================================
      !
      !        this procedure supplies the coefficients a(j), b(j) of the
      !        recurrence relation
      !
      !             b p (x) = (x - a ) p   (x) - b   p   (x)
      !              j j            j   j-1       j-1 j-2
      !
      !        for the various classical (normalized) orthogonal polynomials,
      !        and the 0.0_dp-th moment
      !
      !             muzero = integral w(x) dx
      !
      !        of the given polynomial's weight function w(x).  since the
      !        polynomials are orthonormalized, the tridiagonal matrix is
      !        guaranteed to be symmetric.
      !  *******************************************************************
      !
      !
      !

      !

      !
      integer, intent(in) :: kind
      integer, intent(in) :: n
      integer :: i, nm1
      !
      real (kind = dp), intent(out) :: a(n), b(n), muzero
      real (kind = dp) :: abi
      !
      nm1 = n - 1
      a(1:n) = 0.0_dp
      if (kind.eq.1) then
         !     ================================================================
         !     kind = 1:  legendre polynomials p(x) on (-1, +1), w(x) = 1.
         !     ================================================================
         muzero = 2.0_dp
         do i = 1, nm1
            abi = i
            b(i) = abi / sqrt(4.0_dp * abi * abi - 1.0_dp)
         enddo
      else if (kind.eq.2) then
         !     ================================================================
         !     kind = 2:  chebyshev polynomials of the first kind t(x)
         !                on (-1, +1), w(x) = 1 / sqrt(1 - x*x)
         !     ================================================================
         muzero = pi
         do i = 1, nm1
            b(i) = 0.5_dp
         enddo
         b(1) = sqrt(0.5_dp)
      else if (kind.eq.3) then
         !     ================================================================
         !     kind = 3:  chebyshev polynomials of the second kind u(x)
         !                on (-1, +1), w(x) = sqrt(1 - x*x)
         !     ================================================================
         muzero = pi * 0.5_dp
         do i = 1, nm1
            b(i) = 0.5_dp
         enddo
      else if (kind.eq.4) then
         !     ================================================================
         !     kind = 4:  hermite polynomials h(x) on (-infinity, +infinity),
         !                w(x) = exp(-x**2)
         !     ================================================================
         muzero = sqrt(pi)
         do i = 1, nm1
            b(i) = sqrt(i * 0.5_dp)
         enddo
      else
         print *, 'CLASS', 'kind < 1 or kind > 4', kind
         stop
      endif
      !
   end subroutine class
   !
   !
   subroutine gausq2(n, d, e, z, ierr)
      !  ===================================================================
      !
      !     this subroutine is a translation of an algol procedure,
      !     num. math. 12, 377-383(1968) by martin and wilkinson,
      !     as modified in num. math. 15, 450(1970) by dubrulle.
      !     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
      !     this is a modified version of the 'eispack' routine imtql2.
      !
      !     this subroutine finds the eigenvalues and first comp1.0_dpnts of the
      !     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
      !     method.
      !
      !     on input:
      !
      !        n is the order of the matrix;
      !
      !        d contains the diagonal elements of the input matrix;
      !
      !        e contains the subdiagonal elements of the input matrix
      !          in its first n-1 positions.  e(n) is arbitrary;
      !
      !        z contains the first row of the identity matrix.
      !
      !      on output:
      !
      !        d contains the eigenvalues in ascending order.  if an
      !          error exit is made, the eigenvalues are correct but
      !          unordered for indices 1, 2, ..., ierr-1;
      !
      !        e has been destroyed;
      !
      !        z contains the first comp1.0_dpnts of the orthonormal eigenvectors
      !          of the symmetric tridiagonal matrix.  if an error exit is
      !          made, z contains the eigenvectors associated with the stored
      !          eigenvalues;
      !
      !        ierr is set to
      !          0.0_dp       for normal return,
      !          j          if the j-th eigenvalue has not been
      !                     determined after 30 iterations.
      !
      !  *********************************************************************
      !

      !

      !

      !
      integer, intent(in) :: n
      integer, intent(out) :: ierr
      integer :: i, j, k, l, m, ii, mml
      !
      real (kind = dp), intent(inout) :: d(n), e(n), z(n)
      real (kind = dp) :: b, c, f, g, p, r, s
      real (kind = dp), parameter :: machep = 2.0_dp**(-127)
      !
      !  .....................................................................
      !  Causion: unproper assigned machine parameter could cause
      !           unconvergent problem.
      !  .....................................................................
      !
      ierr = 0
      if (n .eq. 1) return
      !
      e(n) = 0.0_dp
      do l = 1, n
         j = 0
         !     :::::::::: look for small sub-diagonal element ::::::::::
         do
            do m = l, n
               if (m.eq.n) then
                  exit
               else if(abs(e(m)).le.machep * (abs(d(m)) + abs(d(m + 1)))) then
                  exit
               endif
            enddo
            !
            p = d(l)
            if (m .eq. l) then
               exit
            endif
            if (j .eq. 30) then
               !           :::: set error -- no convergence to an
               !                eigenvalue after 30 iterations ::::::::::
               ierr = l
               return
            endif
            j = j + 1
            !        :::::::::: form shift ::::::::::
            g = (d(l + 1) - p) / (2.0_dp * e(l))
            r = sqrt(g * g + 1.0_dp)
            g = d(m) - p + e(l) / (g + sign(r, g))
            s = 1.0_dp
            c = 1.0_dp
            p = 0.0_dp
            mml = m - l
            !
            !        :::::::::: for i=m-1 step -1 until l do -- ::::::::::
            do ii = 1, mml
               i = m - ii
               f = s * e(i)
               b = c * e(i)
               if (abs(f) .ge. abs(g)) then
                  c = g / f
                  r = sqrt(c * c + 1.0_dp)
                  e(i + 1) = f * r
                  s = 1.0_dp / r
                  c = c * s
               else
                  s = f / g
                  r = sqrt(s * s + 1.0_dp)
                  e(i + 1) = g * r
                  c = 1.0_dp / r
                  s = s * c
               endif
               g = d(i + 1) - p
               r = (d(i) - g) * s + 2.0_dp * c * b
               p = s * r
               d(i + 1) = g + p
               g = c * r - b
               !           :::::::::: form first comp1.0_dpnt of vector ::::::::::
               f = z(i + 1)
               z(i + 1) = s * z(i) + c * f
               z(i) = c * z(i) - s * f
            enddo
            !
            d(l) = d(l) - p
            e(l) = g
            e(m) = 0.0_dp
         enddo
      enddo
      !
      !  :::::::::: order eigenvalues and eigenvectors ::::::::::
      do ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
         do j = ii, n
            if (d(j) .lt. p) then
               k = j
               p = d(j)
            endif
         enddo
         if (k .ne. i) then
            d(k) = d(i)
            d(i) = p
            p = z(i)
            z(i) = z(k)
            z(k) = p
         endif
      enddo
      !
   end subroutine gausq2

end module gaussq_mod
