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
  XCBase(int n_spin, int xc_functional[3]);

  XCBase(int n_spin, std::vector<int> xc_functional);

  virtual void evaluate(const std::vector<Real> &r_mesh, double h,
                        const Matrix<Real> &rho_in, int jmt,
                        Matrix<Real> &xc_energy_out, Matrix<Real> &xc_pot_out) = 0;

  virtual void evaluate(const Real rho_in[2], Real &xc_energy_out,
                        Real xc_pot_out[2]) = 0;

  virtual std::string get_name() = 0;
};

}  // namespace lsms

#endif  // LSMS_XCBASE_HPP
