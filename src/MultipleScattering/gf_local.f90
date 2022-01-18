!
!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine gf_local(mtasa, zj_flag, lmax, kkrsz, &
   rins, r_sph, r_mesh, jmt, jws, &
   pnrel, tau00, matom_left, matom_right, &
   zlr_left, zlr_right, jlr, nprpts, &
   ngaussr, &
   cgnt, lmax_cg, &
   dos, dosck, green, dipole, &
   greenIntLLp, &
   ncrit, grwylm, gwwylm, wylm, &
   pi, iprint, istop)
   !     ================================================================
   !
   !
   !     ****************************************************************
   !     input:
   !                tau00
   !                kkrsz   (size of KKR-matrix)
   !                istop   (index of subroutine prog. stops in)
   !     output:
   !                dos     (wigner-seitz cell density of states)
   !                dosck   (muffin-tin density of states)
   !                green   (Green's function)
   !
   !                They are all calculated in the Local frame.
   !     ****************************************************************
   !
   implicit none
   !
   !     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !      include    'atom_param.h'
   !     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   !
   character  istop*32
   character  sname*32
   !
   integer    mtasa
   integer    zj_flag
   integer    lmax
   integer    kkrsz
   integer    jmt
   integer    jws
   integer    m
   integer    ngaussr
   integer    iprint, iprint_dos
   integer lmax_cg, nprpts, ncrit

   integer j, lm1, lm2
   !
   real*8     rins, r_sph, r_mesh(nprpts)
   real*8     cgnt(lmax_cg + 1, (lmax_cg + 1)**2, (lmax_cg + 1)**2)
   real*8     pi
   real*8     sqr2
   real*8 grwylm(*), gwwylm(*)
   !
   complex*16 pnrel
   complex*16 matom_left(lmax + 1)
   complex*16 matom_right(lmax + 1)
   complex*16 tau00(kkrsz, kkrsz)
   complex*16 pzz(kkrsz, kkrsz)
   complex*16 pzj
   complex*16 pzj_full(kkrsz)
   complex*16 pzzck(kkrsz, kkrsz)
   complex*16 pzjck
   complex*16 pzjck_full(kkrsz)
   complex*16 dzz(kkrsz * kkrsz * 3)
   complex*16 dzj
   complex*16 vzz(kkrsz * kkrsz * 3)
   complex*16 vzj
   complex*16 zlr_left(nprpts, 0:lmax)
   complex*16 zlr_right(nprpts, 0:lmax)
   complex*16 jlr(nprpts, 0:lmax)
   complex*16 dos
   complex*16 dosck
   complex*16 green(jws)
   ! dipole(m,1) is the density moment
   ! dipole(m,2) is the gradient of the dipole potential at the origin
   complex*16 dipole(-1:1, 2)
   complex*16 wylm(*)
   !
   complex*16 greenIntLLp(kkrsz, kkrsz)

   parameter (sname = 'gf_local')
   !
   sqr2 = sqrt(2.d0)
   !     =================================================================
   !     call int_zz_zj to calculate:
   !
   !                      Rmt     2            -> ^  ->
   !     pzzck = - 1/pi * int dr*r * 4*pi * Z (r)*Z (r)
   !                       0                 L     L'
   !
   !                      Rmt     2            -> ^  ->
   !     pzjck = + 1/pi * int dr*r * 4*pi * Z (r)*J (r) * delta
   !                       0                 L     L'          LL'
   !
   !                      Rs   3->     -> ^  ->       ->
   !       pzz = - 1/pi * int d r * Z (r)*Z (r)*Sigma(r)
   !                       0         L     L'
   !
   !                      Rs   3->     -> ^  ->       ->
   !       pzj = + 1/pi * int d r * Z (r)*J (r)*Sigma(r) * delta
   !                       0         L     L'                   LL'
   !     =================================================================
   !     -----------------------------------------------------------------
   iprint_dos = iprint

   call int_zz_zj(mtasa, zj_flag, lmax, kkrsz, &
      pnrel, matom_left, matom_right, &
      rins, r_sph, r_mesh, jmt, &
      zlr_left, zlr_right, jlr, nprpts, &
      ngaussr, &
      cgnt, lmax_cg, &
      pzzck, pzz, pzjck, pzj, dzz, dzj, vzz, vzj, &
      pzjck_full, pzj_full, &
      ncrit, grwylm, gwwylm, wylm, &
      iprint_dos, istop)
   !     -----------------------------------------------------------------
   !
   !      write(*,*) "pzz and pzzck after int_zz_zj"
   !      do lm1=1,kkrsz
   !         do lm2=1,kkrsz
   !            write(*,*) lm1, lm2, pzz(lm1,lm2), pzzck(lm1,lm2)
   !         end do
   !      end do
   !      call fstop("gf_local")

   !     calculate radially L, L' resolved Green's function integrated
   !     integrated over the atomic volume

   do lm1 = 1, kkrsz
      do lm2 = 1, kkrsz
         greenIntLLp(lm1, lm2) = tau00(lm1, lm2) * pzz(lm1, lm2)
      end do
   end do
   if (zj_flag .eq. 1) then
      do lm1 = 1, kkrsz
         greenIntLLp(lm1, lm1) = greenIntLLp(lm1, lm1) - pzj_full(lm1)
      end do
   end if


   !     ================================================================
   !     dos =>  ZZ(e)*tau(e) - ZJ(e)....................................
   !     ----------------------------------------------------------------
   call mdosms(zj_flag, kkrsz * kkrsz, dosck, pzzck, pzjck, tau00, pi, &
      iprint, istop)
   !     ----------------------------------------------------------------
   if(iprint.ge.1 .and. zj_flag.eq.1) then
      write(6, '('' e,n(e)_mt:'',2d16.8,2f18.13)') pnrel * pnrel, dosck
   endif
   !     ----------------------------------------------------------------
   call mdosms(zj_flag, kkrsz * kkrsz, dos, pzz, pzj, tau00, pi, &
      iprint_dos, istop)
   !     ----------------------------------------------------------------
   if(iprint.ge.1 .and. zj_flag.eq.1) then
      write(6, '('' e,n(e)_ws:'',2d16.8,2f18.13)') pnrel * pnrel, dos
   endif
   !     ----------------------------------------------------------------
   do m = -1, 1
      call mdosms(zj_flag, kkrsz * kkrsz, dipole(m, 1), &
         dzz(kkrsz * kkrsz * (m + 1) + 1), dzj, tau00, pi, iprint, istop)
      call mdosms(zj_flag, kkrsz * kkrsz, dipole(m, 2), &
         vzz(kkrsz * kkrsz * (m + 1) + 1), vzj, tau00, pi, iprint, istop)
   enddo
   ! change to real harmonic basis so we can take the imag part
   !  dipole is Ylm^* which is
   !        (x+iy)/r2
   !        (  z  )
   !        (-x+iy)/r2
   ! we want to convert it into real harmonics as
   !        (y)
   !        (z)
   !        (x)
   do m = 1, 2
      dzj = (dipole(-1, m) + dipole(1, m)) * dcmplx(0.d0, -1.d0 / sqr2)
      dipole(1, m) = (dipole(-1, m) - dipole(1, m)) / sqr2
      dipole(-1, m) = dzj
   enddo
   !     ----------------------------------------------------------------
   !
   !     ================================================================
   !     green =>  ZZ(r,e)*tau(e) - ZJ(r,e)....................................
   !     ----------------------------------------------------------------
   call mgreen(zj_flag, kkrsz, lmax, jws, nprpts, &
      zlr_left, zlr_right, jlr, green, tau00, &
      iprint, istop)
   !     if(iprint.ge.0)write(6,'(''green'',6d14.6)')green(501),green(jws)
   !     ----------------------------------------------------------------
   !
   !     ================================================================
   if(istop.eq.sname) then
      call fstop(sname)
   else
      return
   endif
   !
end
!
