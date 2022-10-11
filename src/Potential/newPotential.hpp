//
// Created by F.Moitzi on 24.10.2022.
//

#ifndef LSMS_NEWPOTENTIAL_HPP
#define LSMS_NEWPOTENTIAL_HPP

#include "LSMSCommunication.hpp"
#include "Matrix.hpp"
#include "Real.hpp"
#include "SystemParameters.hpp"

extern "C" {

void newexchg_(int *n_spin_pola, Real *sp, Real *rhoup, Real *rhodn, Real *vx,
               Real *enxc, Real *vxout, Real *excout, Real *ro3, Real *dz,
               Real *r_mesh, int *jmt, int *iexch);
}

namespace lsms {

static int counter = 2;

void calculateVRMS(LSMSSystemParameters &lsms, LocalTypeInfo &local);

class Potential {
 public:
  virtual void calculatePotential(LSMSCommunication &comm,
                                  LSMSSystemParameters &lsms,
                                  LocalTypeInfo &local,
                                  CrystalParameters &crystal,
                                  std::vector<Real> &qsub) = 0;
};

class ASAPotential : public Potential {
 public:
  void calculatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms,
                          LocalTypeInfo &local, CrystalParameters &crystal,
                          std::vector<Real> &qsub) override;
};

}  // namespace lsms

#endif  // LSMS_NEWPOTENTIAL_HPP
