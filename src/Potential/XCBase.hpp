//
// Created by F.Moitzi on 24.04.2022.
//

#ifndef LSMS_XCBASE_HPP
#define LSMS_XCBASE_HPP

#include <vector>

#include "Matrix.hpp"
#include "Real.hpp"

namespace lsms {

class XCBase {
 protected:
  int _nSpin;
  std::vector<int> _xcFunctional;

 public:
  XCBase(int nSpin, int xcFunctional[3]);

  XCBase(int nSpin, std::vector<int> xcFunctional);

  virtual void evaluate(std::vector<Real> &rMesh,
                        std::vector<Real> &drMesh,
                        const Matrix<Real> &rhoIn, int jmt,
                        Matrix<Real> &xcEnergyOut, Matrix<Real> &xcPotOut) = 0;

  virtual void evaluate(const Real rhoIn[2], Real &xcEnergyOut,
                        Real xcPotOut[2]) = 0;
};

}  // namespace lsms

#endif  // LSMS_XCBASE_HPP
