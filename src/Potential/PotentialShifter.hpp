#ifndef LSMS_SHIFT_POTENTIALS
#define LSMS_SHIFT_POTENTIALS

#include <vector>

#include "Main/SystemParameters.hpp"
#include "Matrix.hpp"
#include "Real.hpp"

class PotentialShifter {
 public:
  bool vSpinShiftFlag{false};
  double minShift{0.0};
  double maxShift{0.0};

  void resetPotentials(LocalTypeInfo &local);

  void resize(int n);

  void applyShifts(LocalTypeInfo &local);

  void restorePotentials(LocalTypeInfo &local);

  std::vector<Matrix<Real> > vr0{};  // the unshifted potential
};

#endif
