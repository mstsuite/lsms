module xc_mod
   use ISO_C_BINDING
   implicit none

contains

   !
   ! Calculates XC LDA density and potential from the charge density "n".
   !
   subroutine getvxc_scalar(n, relat, c_light, exc, Vxc) bind(C)
      real(c_double), intent(in) :: n ! charge density (scalar)
      real(c_double), intent(in) :: c_light ! speed of light
      logical(c_bool), intent(in) :: relat ! if .true. returns RLDA, otherwise LDA
      real(c_double), intent(out) :: exc ! XC density
      real(c_double), intent(out) :: Vxc ! XC potential

      real(c_double), parameter :: y0 = -0.10498_8
      real(c_double), parameter :: b = 3.72744_8
      real(c_double), parameter :: c = 12.9352_8
      real(c_double), parameter :: A = 0.0621814_8
      real(c_double), parameter :: pi = 3.1415926535897932384626433832795_8

      !     c = 137.0359895

      real(8) :: Q, rs, y, ec, ex, Vc, Vx, beta, mu, R, S

      if (n < epsilon(1.0_8)) then
         exc = 0
         Vxc = 0
         return
      end if

      Q = sqrt(4 * c - b**2)
      rs = (3 / (4 * pi * n))**(1.0_8 / 3)
      y = sqrt(rs)
      ec = A / 2 * (log(y**2 / get_Y(y, b, c)) + 2 * b / Q * atan(Q / (2 * y + b))  &
         - b * y0 / get_Y(y0, b, c) * (&
            log((y - y0)**2 / get_Y(y, b, c)) &
               + 2 * (b + 2 * y0) / Q * atan(Q / (2 * y + b)) &
            ))
      Vc = ec - A / 6 * (c * (y - y0) - b * y0 * y) / ((y - y0) * get_Y(y, b, c))
      ex = -3 / (4 * pi) * (3 * pi**2 * n)**(1.0_8 / 3)
      Vx = 4 * ex / 3

      if (relat) then
         beta = -4 * pi * ex / (3 * c_light)
         mu = sqrt(1 + beta**2)
         R = 1 - 3 * ((beta * mu - log(beta + mu)) / (beta ** 2))**2 / 2
         S = 3 * log(beta + mu) / (2 * beta * mu) - 1.0_8 / 2

         ex = ex * R
         Vx = Vx * S
      end if
      exc = ex + ec
      Vxc = Vx + Vc
   end subroutine

   function get_Y(y, b, c)
      real(8), intent(in) :: y, b, c
      real(8) :: get_Y
      get_Y = y**2 + b * y + c
   end function

end module xc_mod
