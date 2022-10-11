//
// Created by F.Moitzi on 17.01.2023.
//

#include "XCLDA.hpp"

#include "xc_lda_potential.hpp"

namespace lsms {
XCLDA::XCLDA(int n_spin, int *xc_functional) : XCBase(n_spin, xc_functional) {

  rel = true;

}

XCLDA::XCLDA(int n_spin, std::vector<int> xc_functional) : XCBase(n_spin, xc_functional) {
  rel = true;

}

void XCLDA::evaluate(const Real *rho_in, Real &xc_energy_out, Real *xc_pot_out) {

}

std::string XCLDA::get_name() {
  return "LDA";
}

void XCLDA::evaluate(const std::vector<Real> &r_mesh,
                     double h,
                     const Matrix<Real> &rho_in,
                     int jmt,
                     Matrix<Real> &xc_energy_out,
                     Matrix<Real> &xc_pot_out) {

  if (rel) {

    for (int i = 0; i < jmt; i++) {

      VxcRLDA(r_mesh[i],
              rho_in(i, 0) / (4.0 * r_mesh[i] * r_mesh[i] * M_PI),
              xc_pot_out(i, 0),
              xc_energy_out(i, 0)
      );

    }

  } else {

    for (int i = 0; i < jmt; i++) {

      VxcLDA(r_mesh[i],
             rho_in(i, 0) / (4.0 * r_mesh[i] * r_mesh[i] * M_PI),
             xc_pot_out(i, 0),
             xc_energy_out(i, 0)
      );

    }

  }

}

} // lsms