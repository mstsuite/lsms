//
// Created by F.Moitzi on 17.01.2023.
//

#ifndef LSMS_SRC_POTENTIAL_XCLDA_HPP_
#define LSMS_SRC_POTENTIAL_XCLDA_HPP_

#include "XCBase.hpp"

#include <vector>

namespace lsms {

class XCLDA : public XCBase {

 private:

  bool rel;

 public:
  XCLDA(int n_spin, int xc_functional[3]);

  XCLDA(int n_spin, std::vector<int> xc_functional);

  void evaluate(const std::vector<Real> &r_mesh, double h,
                const Matrix<Real> &rho_in, int jmt,
                Matrix<Real> &xc_energy_out, Matrix<Real> &xc_pot_out) override;

  void evaluate(const Real rho_in[2], Real &xc_energy_out,
                Real xc_pot_out[2]) override;

  std::string get_name() override;

};

} // lsms

#endif //LSMS_SRC_POTENTIAL_XCLDA_HPP_
