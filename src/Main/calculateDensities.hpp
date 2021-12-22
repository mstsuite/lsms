#ifndef LSMS_CALCULATEDENSITIES_HPP
#define LSMS_CALCULATEDENSITIES_HPP

#include <vector>

#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"
#include "Array3d.hpp"

#include "SystemParameters.hpp"
#include "checkConsistency.hpp"

void calculateChargeDensity(LSMSSystemParameters &lsms, AtomData &atom, Real edote,
                            Matrix<Real> &rhonew, Real &qvalmt, Real *qrms);

void calculateAllLocalChargeDensities(LSMSSystemParameters &lsms, LocalTypeInfo &local);

void calculateDensities(LSMSSystemParameters &lsms, int iatom, int is, int ie, int nume, Complex energy, Complex dele1,
                        Matrix<Complex> &dos, Matrix<Complex> &dosck, Array3d<Complex> &green,
                        Array3d<Complex> &dipole,
                        AtomData &atom);

void checkAllLocalCharges(LSMSSystemParameters &lsms, LocalTypeInfo &local);

void calculateLocalQrms(LSMSSystemParameters &lsms, LocalTypeInfo &local);


extern "C"
{
  void interp_(Real *r, Real *f, int *nr, Real *rs, Real *ps, Real *dps, int *deriv);
  void newint_(int *nr, Real *r, Real *f, Real *g, int *ip0);
}

#endif
