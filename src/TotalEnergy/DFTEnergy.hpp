//
// Created by F.Moitzi on 05.07.2022.
//

#ifndef LSMS_DFTENERGY_HPP
#define LSMS_DFTENERGY_HPP

#include "LSMSCommunication.hpp"
#include "Real.hpp"

namespace lsms {

class DFTEnergy {
 public:
  DFTEnergy() = default;

  Real zero_point = 0.0;

  Real core_eigen = 0.0;

  Real semicore_eigen = 0.0;

  Real one_ele = 0.0;

  Real ks = 0.0;

  Real kinetic = 0.0;

  Real hartree = 0.0;

  Real core_hartree = 0.0;

  Real coloumb = 0.0;

  Real xc = 0.0;

  Real lsf = 0.0;

  Real total = 0.0;

  Real madelung = 0.0;

  Real it_madelung = 0.0;

  Real it_xc = 0.0;

  Real mtz = 0.0;

  Real u0 = 0.0;

  DFTEnergy &operator+=(const DFTEnergy &rhs) {
    this->zero_point += rhs.zero_point;
    this->core_eigen += rhs.core_eigen;
    this->semicore_eigen += rhs.semicore_eigen;
    this->one_ele += rhs.one_ele;
    this->ks += rhs.ks;
    this->kinetic += rhs.kinetic;
    this->hartree += rhs.hartree;
    this->core_hartree += rhs.core_hartree;
    this->coloumb += rhs.coloumb;
    this->xc += rhs.xc;
    this->lsf += rhs.lsf;
    this->total += rhs.total;

    /**
     * There are global quantities
     */

    return *this;
  }
};

void globalSum(LSMSCommunication &comm, DFTEnergy &dft_energy);

void print_dft_energy(const DFTEnergy &energy);

}  // namespace lsms

#endif  // LSMS_DFTENERGY_HPP
