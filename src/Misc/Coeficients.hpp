#ifndef LSMS_COEFICIENTS_HPP
#define LSMS_COEFICIENTS_HPP

#include <cmath>
#include <vector>

#include "Array3d.hpp"
#include "Complex.hpp"
#include "Indices.hpp"
#include "Main/SystemParameters.hpp"
#include "Matrix.hpp"

extern "C" {
void cgaunt_(int *lmax, Real *clm, Real *plmg, Real *tg, Real *wg, Real *cgnt,
             int *lofk, int *mofk, int *iprint, char *istop);

void ifacts_(int *lmax, Complex *illp, Complex *ilp1, int *iprint, char *istop,
             int istop_len);
}

class SphericalHarmonicsCoeficients {
 public:
  static int lmax;
  static std::vector<Real> clm;

  static void init(int _lmax) {
    lmax = _lmax;
    clm.resize((lmax + 1) * (lmax + 2) / 2);
    // coefficients for normalized associated Legendre functions (P_00 =
    // \sqrt{1/4pi})
    for (int lm = 0; lm < (lmax + 1) * (lmax + 2) / 2; lm++) {
      clm[lm] = 1.0;
    }
  }
};

class GauntCoeficients {
 public:
  static int lmax;
  static Array3d<Real> cgnt;

  static void init(LSMSSystemParameters &lsms) {
    lmax = AngularMomentumIndices::lmax / 2;
    // lmax=lsms.maxlmax;

    //    if(useNewGaunt)
    //    {
    //      calculateGauntCoeficients(lmax, cgnt, a);
    //      if(testNewGaunt)
    //      {
    //        const Real tol=1.0e-12;
    //        std::vector<Real> tg,wg;
    //        Matrix<Real> plmg;
    //        Array3d<Real> cgntTest;
    //        tg.resize(2*(2*lmax+1)); wg.resize(2*(2*lmax+1));
    //        cgntTest.resize(lmax+1,(lmax+1)*(lmax+1),(lmax+1)*(lmax+1));
    //        plmg.resize(((2*lmax+1)*(2*lmax+2))/2,2*lmax+1);
    //
    //        cgaunt_(&lmax,&SphericalHarmonicsCoeficients::clm[0],&plmg(0,0),&tg[0],&wg[0],
    //                &cgntTest(0,0,0),
    //                &AngularMomentumIndices::lofk[0],&a.mofk[0],
    //                &lsms.global.iprint,lsms.global.istop);
    //
    //        for(int l1=0; l1<lmax+1; l1++)
    //          for(int j2=0; j2<(lmax+1)*(lmax+1); j2++)
    //            for(int j3=0; j3<(lmax+1)*(lmax+1); j3++)
    //            {
    //              if(std::abs(cgnt(l1,j2,j3)-cgntTest(l1,j2,j3))>tol)
    //              {
    //                printf("Difference in Gaunt Coeficients: %d %d %d :
    //                new=%.12lf old=%.12lf\n",
    //                       l1,j2,j3,cgnt(l1,j2,j3),cgntTest(l1,j2,j3));
    //                exit(1);
    //              }
    //            }
    //      }
    //    } else {

    std::vector<Real> tg, wg;
    Matrix<Real> plmg;
    tg.resize(2 * (2 * lmax + 1));
    wg.resize(2 * (2 * lmax + 1));
    cgnt.resize(lmax + 1, (lmax + 1) * (lmax + 1), (lmax + 1) * (lmax + 1));
    plmg.resize(((2 * lmax + 1) * (2 * lmax + 2)) / 2, 2 * lmax + 1);

    cgaunt_(&lmax, &SphericalHarmonicsCoeficients::clm[0], &plmg(0, 0), &tg[0],
            &wg[0], &cgnt(0, 0, 0), &AngularMomentumIndices::lofk[0],
            &AngularMomentumIndices::mofk[0], &lsms.global.iprint,
            lsms.global.istop);
    // }
  }
};

class IFactors {
 public:
  static int lmax;
  static Matrix<Complex> illp;
  static std::vector<Complex> ilp1;

  static void init(LSMSSystemParameters &lsms, int _lmax) {
    lmax = _lmax;
    ilp1.resize(2 * _lmax + 1);
    illp.resize((_lmax + 1) * (_lmax + 1), (_lmax + 1) * (_lmax + 1));
    ifacts_(&_lmax, &illp(0, 0), &ilp1[0], &lsms.global.iprint,
            lsms.global.istop, 32);
  }
};

#endif
